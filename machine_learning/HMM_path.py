#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-11-06
# Github acct: edevog
import sys
import pandas as pd

'''
This program decodes a Hidden Markov Model to find the most probable hidden path
for a given sequence.

Input: A text file with a sequence, a list of unique characters that comprise
the sequence, a list of hidden states, a probability matrix of transitioning
from one state to another, a probability matrix of being in a particular state
for a particular character in the sequence. All matrixes and lists are tab
delimited.
Ex:
xyxzzxyxyy
--------
x   y   z
--------
A   B
--------
    A   B
A   0.641   0.359
B   0.729   0.271
--------
    x   y   z
A   0.117   0.691   0.192
B   0.097   0.42    0.483

Output:
The most probable list of states for our given sequence and model
AAABBAAAAA
'''

class Path:
    '''
    Calculates and defines values related to a Hidden Markov Model

    Attributes:
        seq = the emitted sequence from the HMM
        states = the states in the HMM
        transition = the data frame of the probabilities of transitioning
                        between states
        emission = the data from of the probabilities of being in a certain
                    state for a certain symbol in the sequence
    '''
    def __init__(self, seq, states, transition, emission):
        self.seq = seq
        self.states = states
        self.transition = transition
        self.emission = emission


    def findPath(self, seq, states, transition, emission):
        '''
        Finds the most probable path of states given an emitted sequence from
        a Hidden Markov model
        '''
        #lists the positions in the sequence to use as references
        pos = list(range(0, len(self.seq)))
        #an empty list to later store initial values of the score matrix
        s_m = list()
        #an empty list to later store initial values of the path matrix
        p_m = list()
        #iterates through the number of states to create an empty matrix of the
        #appropriate size for the score matrix and path matrix
        for i in range(len(states)):
            s_m.append([0.0]*len(seq))
            p_m.append([""]*len(seq))
        #defining the score and path  matrices as data frames so it is more
        #intuitive to call values from the matrices
        score_matrix = pd.DataFrame(s_m, columns = pos, index = states)
        path_matrix = pd.DataFrame(p_m, columns = pos, index = states)
        #iterate through the sequence and calculates the score for
        #each state
        #score = max(state_i-1 * emit(syb_i, state_i)
                #* tran(state_i-1, state_i)) for all states_i-1
        for i in pos:
            #the first symbol in the emitted sequence is treated differently
            #because it has no preivous scores to calculate its value
            if i == 0:
                for s in states:
                    #transition value = 1/(#of states)
                    #prev_score = 1
                    score_matrix[i][s] = (1/len(states))*emission[seq[i]][s]
            else:
                #iterates through the states we could go into next
                for next in states:
                    #iterates through the states we could have come from
                    for prev in states:
                        #defining the values of each part of the
                        #score calculation to make the code more readable

                        #finds the score of one of the previous state
                        prev_score = score_matrix[i-1][prev]
                        tran = transition[next][prev]
                        emit = emission[seq[i]][next]
                        score = prev_score * emit * tran
                        #check if the current score is better than the
                        #previously calculated score for this symbol and state
                        if score > score_matrix[i][next]:
                            #if it's better store the value and find its path
                            score_matrix[i][next] = score
                            path_matrix[i-1][next] = prev
        f_score = 0
        f_path = ""
        #iterate through the states in the final sequence value
        for s in states:
            #store the larger score and its path
            if f_score < score_matrix[len(seq)-1][s]:
                f_score = score_matrix[len(seq)-1][s]
                # f_path = path_matrix[len(seq)-1][s]
                f_path = [s]
                place = int(len(self.seq) - 2)
                while len(f_path) < len(self.seq):
                    f_path.append(path_matrix[place][s])
                    s = path_matrix[place][s]
                    place -= 1
        f_path.reverse()
        return f_path



def main():
    #the first list of data file is the sequence
    seq = list(sys.stdin.readline().strip())
    #the rest of our data file that we will continue to parse
    f = sys.stdin.readlines()
    #the list of unique symbols that comprise our sequence
    sybs = f[1].strip().split()
    #the list of states
    states = f[3].strip().split()
    #reading a transition matrix with dimensions the same length as the
    #number of states in our model
    tran_matrix = list()
    for i in range(len(states)):
        #using pop so the code is dynamic for the number of states
        tran_matrix.append(list(map(float, f.pop(6).split()[1:])))
    #defines the transition matrix as a dataframe with the column and row
    #names as the states so it is more intuitive to call data from the matrix
    transition = pd.DataFrame(tran_matrix, columns = states, index = states)
    #reading in an emission matrix with dimensions the same length as
    #the number of states in our model
    emit_matrix = list()
    for i in range(len(states)):
        #using pop so the code is dynamic for the number of states
        emit_matrix.append(list(map(float, f.pop(8).split()[1:])))
    emission = pd.DataFrame(emit_matrix, columns = sybs, index = states)


    my_path = Path(seq, states, transition, emission)
    #calculates the most probable path of states given the emitted sequence
    f_path = my_path.findPath(seq, states, transition, emission)
    #prints the path
    print("".join(f_path))


if __name__ == '__main__':
    main()
