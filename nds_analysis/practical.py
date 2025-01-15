"""
This module contains the implementation of 'practical' NDS.

The practical NDS is a simple implementation of the NDS algorithm running on the Left Move Problem.

"""

import random
import copy

SAMPLING_SIZE = 100000


def playout(rollout, h):
    while True:
        if len(rollout) == h:
            return rollout
        rollout = rollout + [random.randint(0, 1)]


def practical_nmcs(rollout, level, h, stop=None):
    if stop is None:
        stop = h
    if level == 0:
        return playout(rollout, h)
    else:
        best = []
        while True:
            if len(rollout) == h or len(rollout) == stop:
                return rollout
            left = copy.deepcopy(rollout)
            left = left + [1]
            left = practical_nmcs(left, level - 1, h)
            right = copy.deepcopy(rollout)
            right = right + [0]
            right = practical_nmcs(right, level - 1, h)
            local = left
            if random.randint(0, 1) == 0:
                local = right
            if sum(left) > sum(right):
                local = left
            if sum(right) > sum(left):
                local = right
            if sum(local) > sum(best):
                best = local
            rollout = rollout + [best[len(rollout)]]


def practical_nds_depth_two(rollout, level, h):
    if level == 0:
        return playout(rollout, h)
    else:
        best = []
        while True:
            if len(rollout) == h:
                return rollout
            leftleft = copy.deepcopy(rollout)
            leftleft = leftleft + [1, 1]
            leftleft = practical_nds_depth_two(leftleft, level - 1, h)
            leftright = copy.deepcopy(rollout)
            leftright = leftright + [1, 0]
            leftright = practical_nds_depth_two(leftright, level - 1, h)
            rightleft = copy.deepcopy(rollout)
            rightleft = rightleft + [0, 1]
            rightleft = practical_nds_depth_two(rightleft, level - 1, h)
            rightright = copy.deepcopy(rollout)
            rightright = rightright + [0, 0]
            rightright = practical_nds_depth_two(rightright, level - 1, h)
            local = []
            if (
                sum(leftleft) >= sum(leftright)
                and sum(leftleft) >= sum(rightleft)
                and sum(leftleft) >= sum(rightright)
            ):
                local.append(leftleft)
            if (
                sum(leftright) >= sum(leftleft)
                and sum(leftright) >= sum(rightleft)
                and sum(leftright) >= sum(rightright)
            ):
                local.append(leftright)
            if (
                sum(rightleft) >= sum(leftleft)
                and sum(rightleft) >= sum(rightright)
                and sum(rightleft) >= sum(leftright)
            ):
                local.append(rightleft)
            if (
                sum(rightright) >= sum(leftleft)
                and sum(rightright) >= sum(rightleft)
                and sum(rightright) >= sum(leftright)
            ):
                local.append(rightright)
            i = random.randint(0, len(local) - 1)
            if sum(local[i]) > sum(best):
                best = local[i]
            rollout = rollout + [best[len(rollout)], best[len(rollout) + 1]]


def practical_nds_depth_two_step_one(rollout, level, h):
    rollout = copy.deepcopy(rollout)
    if level == 0:
        return playout(rollout, h)
    else:
        best = []
        while True:
            if len(rollout) == h:
                return rollout
            if len(rollout) + 2 > h:
                left = copy.deepcopy(rollout)
                left = left + [1]
                left = practical_nmcs(left, level - 1, h)
                right = copy.deepcopy(rollout)
                right = right + [0]
                right = practical_nmcs(right, level - 1, h)
                local = left
                if random.randint(0, 1) == 0:
                    local = right
                if sum(left) > sum(right):
                    local = left
                if sum(right) > sum(left):
                    local = right
                if sum(local) > sum(best):
                    best = local
                rollout += [best[len(rollout)]]
            else:
                if len(rollout) == h:
                    return rollout
                leftleft = copy.deepcopy(rollout)
                leftleft = leftleft + [1, 1]
                leftleft = practical_nds_depth_two_step_one(leftleft, level - 1, h)
                leftright = copy.deepcopy(rollout)
                leftright = leftright + [1, 0]
                leftright = practical_nds_depth_two_step_one(leftright, level - 1, h)
                rightleft = copy.deepcopy(rollout)
                rightleft = rightleft + [0, 1]
                rightleft = practical_nds_depth_two_step_one(rightleft, level - 1, h)
                rightright = copy.deepcopy(rollout)
                rightright = rightright + [0, 0]
                rightright = practical_nds_depth_two_step_one(rightright, level - 1, h)
                max_sum = max(
                    [sum(leftleft), sum(leftright), sum(rightleft), sum(rightright)]
                )
                local = [
                    elem
                    for elem in [leftleft, leftright, rightleft, rightright]
                    if sum(elem) == max_sum
                ]
                i = random.randint(0, len(local) - 1)
                if sum(local[i]) > sum(best):
                    best = local[i]
                rollout = rollout + [best[len(rollout)]]


def run_practical_nested(level, depth, step, height):
    func_call = None
    if depth == 1:
        func_call = practical_nmcs
    elif depth == 2:
        if step == 2:
            func_call = practical_nds_depth_two
        elif step == 1:
            func_call = practical_nds_depth_two_step_one
        else:
            print(
                "Practical nested "
                + str(level)
                + ", depth "
                + str(depth)
                + ", step "
                + str(step)
                + ", height "
                + str(height)
                + " => "
                + "Not implemented"
            )
            return
    else:
        print(
            "Practical nested "
            + str(level)
            + ", depth "
            + str(depth)
            + ", step, "
            + str(step)
            + ", height "
            + str(height)
            + " => "
            + "Not implemented"
        )
        return
    score_counter = {}
    for s in range(height + 1):
        score_counter[s] = 0
    average = 0
    for _ in range(SAMPLING_SIZE):
        rollout = func_call([], level, height)
        score = sum(rollout)
        if score in score_counter:
            score_counter[score] += 1
        else:
            score_counter[score] = 1
        average += score
    average = average / SAMPLING_SIZE
    for key in score_counter:
        score_counter[key] = score_counter[key] / SAMPLING_SIZE
    print(
        "Practical nested "
        + str(level)
        + ", depth "
        + str(depth)
        + ", step, "
        + str(step)
        + ", height "
        + str(height)
        + " => "
        + str(average)
    )
    return score_counter
