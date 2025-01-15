"""
This module contains the implementation of 'theoretical' NDS.

The theoretical NDS is an implementation of the NDS algorithm running on the Left Move Problem which computes the
probability of each sequence of moves to be selected by the NDS algorithm.

"""

import numpy as np

DEBUG = False


def generate_sub_sequences(sequence_set, depth):
    previous_sequence_set = [([], sequence_set, [])]
    current_split = 0
    while current_split < depth:
        new_sequence_set = []
        for hist, sequences, _ in previous_sequence_set:
            left_sequences = []
            right_sequences = []
            for playout in range(len(sequences)):
                if sequences[playout][0] == 1:
                    left_sequences.append(sequences[playout][1:])
                else:
                    right_sequences.append(sequences[playout][1:])
            new_sequence_set.append(
                (hist + [1], left_sequences, [0.0] * len(left_sequences))
            )
            new_sequence_set.append(
                (hist + [0], right_sequences, [0.0] * len(right_sequences))
            )
        previous_sequence_set = list(new_sequence_set)
        current_split += 1
    return new_sequence_set


def proba_sampling_theory(sequence_set, level, depth, step):
    if DEBUG:
        print(
            "|\t" * nesting,
            "proba_sampling_theory",
            "level",
            level,
            "depth",
            depth,
            "step",
            step,
            "sequence_set",
            sequence_set,
        )
    sequence_set = [list(sequence) for sequence in sequence_set]
    depth = min(depth, len(sequence_set[0]))

    if level >= len(sequence_set[0]):
        return [1.0] + [0.0] * (len(sequence_set) - 1)

    new_sequence_set = generate_sub_sequences(sequence_set, depth)
    if level > 1:
        results = []
        for hist, sequences, _ in new_sequence_set:
            results.append(
                (
                    hist,
                    sequences,
                    nested_depth_theory(level - 1, depth, step, sequences, hist, [])[0][
                        2
                    ],
                )
            )
        sequence_list_iterator_list = [
            iter(sequence_list) for _, sequence_list, _ in results
        ]
        add_score = [sum(hist) for hist, _, _ in results]
        old_probabilities = [probabilities for _, _, probabilities in results]
        new_probabilities = [probabilities for _, _, probabilities in new_sequence_set]
        sequence_list = [
            next(sequence_list) for sequence_list in sequence_list_iterator_list
        ]
        done = False
        idx_sequence_list = [0] * len(sequence_list)
        while not done:
            current_scores = [sum(seq) for seq in sequence_list]
            current_scores = [
                score + add_score[i] for i, score in enumerate(current_scores)
            ]

            max_score = max(current_scores)
            max_score_idxs = [
                i for i, score in enumerate(current_scores) if score == max_score
            ]
            old_proba = [
                old_probabilities[i][idx_sequence_list[i]]
                for i, _ in enumerate(current_scores)
            ]
            for i, score in enumerate(current_scores):
                if score == max_score:
                    new_probabilities[i][idx_sequence_list[i]] += np.prod(
                        old_proba
                    ) / len(max_score_idxs)

            for i in range(len(sequence_list_iterator_list) - 1, -1, -1):
                try:
                    sequence_list[i] = next(sequence_list_iterator_list[i])
                    idx_sequence_list[i] += 1
                    break
                except StopIteration:
                    if i == 0:
                        done = True
                        break
                    sequence_list_iterator_list[i] = iter(new_sequence_set[i][1])
                    sequence_list[i] = next(sequence_list_iterator_list[i])
                    idx_sequence_list[i] = 0

        result = [[], []]
        for hist, sequences, probabilities in new_sequence_set:
            for sequence, proba in zip(sequences, probabilities):
                result[0].append(hist + sequence)
                result[1].append(proba)

        return result[1]
    else:
        probabilities = [probabilities for _, _, probabilities in new_sequence_set]
        sums = [
            [sum(sequence) for sequence in sequence_list]
            for _, sequence_list, _ in new_sequence_set
        ]

        max_score = len(sequence_set[0])
        sequence_list_scores_counter_pointer_list = []
        for i in range(len(new_sequence_set)):

            hist_score = sum(new_sequence_set[i][0])
            score_counter_pointer = []
            for _ in range(max_score + 1):
                score_counter_pointer.append([])

            for j in range(len(new_sequence_set[i][1])):
                score_counter_pointer[sums[i][j] + hist_score].append(j)
            sequence_list_scores_counter_pointer_list.append(score_counter_pointer)

        for current_set, sequence_list_scores_counter_pointer in enumerate(
            sequence_list_scores_counter_pointer_list
        ):
            for current_score, score_counter_pointer in enumerate(
                sequence_list_scores_counter_pointer
            ):
                if len(score_counter_pointer) == 0:
                    continue

                probability_update = 0
                combination = []

                impossible_combination = False
                for other_set, scp in enumerate(
                    sequence_list_scores_counter_pointer_list
                ):
                    if other_set == current_set:
                        continue

                    subset = [
                        (sc_score, len(sc))
                        for sc_score, sc in enumerate(scp)
                        if sc_score <= current_score and len(sc) > 0
                    ]
                    if len(subset) == 0:
                        impossible_combination = True
                        break
                    combination.append(subset)

                if impossible_combination:
                    continue

                combination_iterator = [iter(c) for c in combination]
                combination_values = [next(c_itr) for c_itr in combination_iterator]

                done = False
                while not done:

                    nb_current_score = []
                    nb_smaller_score = []
                    for score, len_set in combination_values:
                        if score == current_score:
                            nb_current_score.append(len_set)
                        else:
                            nb_smaller_score.append(len_set)

                    if len(nb_current_score) == 0:
                        probability_update += np.prod(nb_smaller_score)
                    else:
                        probability_update += (
                            np.prod(nb_smaller_score)
                            * np.prod(nb_current_score)
                            * (1 / (1 + len(nb_current_score)))
                        )

                    for i in range(len(combination_values) - 1, -1, -1):
                        try:
                            combination_values[i] = next(combination_iterator[i])
                            break
                        except StopIteration:
                            if i == 0:
                                done = True
                                break
                            combination_iterator[i] = iter(combination[i])
                            combination_values[i] = next(combination_iterator[i])

                for sequence_idx in score_counter_pointer:
                    probabilities[current_set][sequence_idx] += probability_update

        for _, _, probabilities in new_sequence_set:
            for i in range(len(probabilities)):
                probabilities[i] /= (len(sequence_set) / (2**depth)) ** (2**depth)

        result = [[], []]
        for hist, sequences, probabilities in new_sequence_set:
            for sequence, proba in zip(sequences, probabilities):
                result[0].append(hist + sequence)
                result[1].append(proba)

        return result[1]


def proba_selection_theory(sequence_set, sampling_probabilities, proba_dist):
    probabilites = [0.0] * len(sequence_set)
    scores = [sum(sequence_set[i]) for i in range(len(sequence_set))]
    for i in range(len(sampling_probabilities)):
        for j in range(len(proba_dist)):
            if scores[i] > scores[j]:
                probabilites[i] += (
                    sampling_probabilities[i] * proba_dist[j]
                    + sampling_probabilities[j] * proba_dist[i]
                )
            elif scores[i] == scores[j]:
                probabilites[i] += sampling_probabilities[j] * proba_dist[i]
    return probabilites


def proba_combine_theory(sequence_set, proba_dist, level, depth, step):
    sampling_probabilities = proba_sampling_theory(sequence_set, level, depth, step)
    combined_probabilities = proba_selection_theory(
        sequence_set, sampling_probabilities, proba_dist
    )
    return combined_probabilities


def split_sequence_depth_theory(sequence_set_proba_dist_pair):
    new_history_sequence_set_proba_dist_pair = []
    for history, sequence_set, proba_dist in sequence_set_proba_dist_pair:
        left_sequences = []
        left_probabilities = []
        right_sequences = []
        right_probabilities = []
        for playout in range(len(sequence_set)):
            if sequence_set[playout][0] == 1:
                left_sequences.append(sequence_set[playout][1:])
                left_probabilities.append(proba_dist[playout])
            else:
                right_sequences.append(sequence_set[playout][1:])
                right_probabilities.append(proba_dist[playout])
        new_history_sequence_set_proba_dist_pair.append(
            (history + [1], left_sequences, left_probabilities)
        )
        new_history_sequence_set_proba_dist_pair.append(
            (history + [0], right_sequences, right_probabilities)
        )
    return new_history_sequence_set_proba_dist_pair


def merge_sequence_depth_theory(sequence_set_proba_dist_pair):
    new_history_sequence_set_proba_dist_pair = []
    for i in range(0, len(sequence_set_proba_dist_pair), 2):
        left_history, left_sequence_set, left_proba_dist = sequence_set_proba_dist_pair[
            i
        ]
        right_history, right_sequence_set, right_proba_dist = (
            sequence_set_proba_dist_pair[i + 1]
        )
        assert left_history == right_history
        assert len(left_sequence_set) == len(right_sequence_set)
        assert left_sequence_set[0][:-1] == right_sequence_set[0][:-1]
        new_history_sequence_set_proba_dist_pair.append(
            (
                left_history,
                [seq[:-1] for seq in left_sequence_set],
                [
                    left_proba_dist[i] + right_proba_dist[i]
                    for i in range(len(left_proba_dist))
                ],
            )
        )
    return new_history_sequence_set_proba_dist_pair


def assemble_sequence_depth_theory(sequence_set_proba_dist_pair):
    new_history_sequence_set_proba_dist_pair = [
        (sequence_set_proba_dist_pair[0][0], [], [])
    ]
    for _, sequence_set, proba_dist in sequence_set_proba_dist_pair:
        for sequence, proba in zip(sequence_set, proba_dist):
            new_history_sequence_set_proba_dist_pair[0][1].append(sequence)
            new_history_sequence_set_proba_dist_pair[0][2].append(proba)
    return new_history_sequence_set_proba_dist_pair


nesting = 0


def nested_depth_theory(level, depth, step, play_set, current_history, prob_mem):
    global nesting
    nesting += 1
    if DEBUG:
        print(
            "|\t" * nesting,
            "nested_depth_theory",
            "level",
            level,
            "depth",
            depth,
            "step",
            step,
            "play_set",
            play_set,
            "current_history",
            current_history,
            "prob_mem",
            prob_mem,
        )

    if len(play_set[0]) == 0:
        nesting -= 1
        return [([], [current_history], prob_mem)]
    if prob_mem == []:
        new_probabilities = proba_sampling_theory(play_set, level, depth, step)
    else:
        new_probabilities = proba_combine_theory(play_set, prob_mem, level, depth, step)

    current_split = 0
    history_sequence_set_proba_dist_tuples = [([], play_set, new_probabilities)]
    while current_split < step:
        history_sequence_set_proba_dist_tuples = split_sequence_depth_theory(
            history_sequence_set_proba_dist_tuples
        )
        current_split += 1

    new_history_sequence_set_proba_dist_tuples = []
    for history, sequence_set, proba_dist in history_sequence_set_proba_dist_tuples:
        new_history_sequence_set_proba_dist_tuples.extend(
            nested_depth_theory(
                level, depth, step, sequence_set, current_history + history, proba_dist
            )
        )

    new_history_sequence_set_proba_dist_tuples = assemble_sequence_depth_theory(
        new_history_sequence_set_proba_dist_tuples
    )
    nesting -= 1
    return new_history_sequence_set_proba_dist_tuples


def sequences(h, s, l):
    if h == 0:
        l.append(s)
        return
    sequences(h - 1, s + [1], l)
    sequences(h - 1, s + [0], l)


def run_theory_nested(level, depth, step, height):
    seq = []
    sequences(height, [], seq)
    out = nested_depth_theory(level, depth, step, seq, [], [])
    _, s, p = out[0]
    score_counter = {}
    average = 0
    for sp, prob in zip(s, p):
        score = sum(sp)
        if score in score_counter:
            score_counter[score] += prob
        else:
            score_counter[score] = prob
        average += sum(sp) * prob
    print(
        "Theory Search nested "
        + str(level)
        + ", depth "
        + str(depth)
        + ", step "
        + str(step)
        + ", height "
        + str(height)
        + " => "
        + str(average)
    )
    return score_counter