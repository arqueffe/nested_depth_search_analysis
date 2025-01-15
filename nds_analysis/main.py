"""
This module contains the implementation of the graphing functions for the NDS analysis.

"""

import matplotlib.pyplot as plt
from practical import *
from theory import *

if __name__ == '__main__':
    
    
    height = 4
    theoretical_score_table = {}
    practical_score_table = {}
    for level in range(1, 3, 1):
        for depth in range(1, 2, 1):
            if height % depth != 0:
                continue
            for step in range(1, 3, 1):
                if  height % step * level != 0: #step > depth or
                    continue
                try:
                    practical_score_table[(level, depth, step)] = run_practical_nested(level, depth, step, height)
                    print(practical_score_table[(level, depth, step)])
                except Exception as e:
                    # print(e)
                    print("Practical nested ", level, depth, step, height, " => Failed")
                try:
                    theoretical_score_table[(level, depth, step)] = run_theory_nested(level, depth, step, height)
                    print(theoretical_score_table[(level, depth, step)])
                except Exception as e:
                    # print(e)
                    print("Theory nested ", level, depth, step, height, " => Failed")
                    
    print(theoretical_score_table)
    print(practical_score_table)
    # Graph side by the side the practical and theoretical results with y the probability is x the score
    for key in theoretical_score_table:
        (level, depth, step) = key
        practical_scores = practical_score_table[key]
        theoretical_scores = theoretical_score_table[key]
        plt.plot(practical_scores.keys(), practical_scores.values(), label="Practical " + str(key))
        plt.plot(theoretical_scores.keys(), theoretical_scores.values(), label="Theoretical " + str(key))
        mean_squared_error = []
        for score in practical_scores:
            mean_squared_error.append((practical_scores[score] - theoretical_scores[score]) ** 2)
        print("Mean Squared Error: ", sum(mean_squared_error) / len(mean_squared_error))
        root_mean_squared_error = []
        for score in practical_scores:
            root_mean_squared_error.append((practical_scores[score] - theoretical_scores[score]) ** 2)
        root_mean_squared_error = (sum(root_mean_squared_error) / len(root_mean_squared_error)) ** 0.5
        print("Root Mean Squared Error: ", root_mean_squared_error)
        mean_absolute_error = []
        for score in practical_scores:
            mean_absolute_error.append(abs(practical_scores[score] - theoretical_scores[score]))
        mean_absolute_error = sum(mean_absolute_error) / len(mean_absolute_error)
        print("Mean Absolute Error: ", mean_absolute_error)
        plt.xlabel("Score")
        plt.ylabel("Probability")
        plt.title("Practical vs Theoretical")
        plt.legend()
    plt.show()
    