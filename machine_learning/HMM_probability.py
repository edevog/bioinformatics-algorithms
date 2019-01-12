#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-11-06
# Github acct: edevog
import sys
import numpy as np
import pandas as pd
'''
This program calculates the probability of a string emitted by a HMM.

Input: A text file with the emitted string, the unique characters from the
string, the possible states, a matrix of probability of trasitioning from one
state to another, and an emission matrix of the probabilities of a certain
character from string being in a certain state.
Ex:
xzyyzzyzyy
--------
x   y   z
--------
A   B
--------
    A   B
A   0.303   0.697
B   0.831   0.169
--------
    x   y   z
A   0.533   0.065   0.402
B   0.342   0.334   0.324

Output:
The probability of the string emitted by the HMM.
Ex:
1.1005510319694847e-06
'''

class HMM:
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


    def calcProb(self):
        '''
        Calculates the probability of emitting a sequence from the HMM based
        on the transition and emission matrices.
        '''
        #lists the positions in the sequence to use as references
        pos = list(range(0, len(self.seq)))
        #an empty list to later store initial values of the score matrix
        s_m = list()
        #an empty list to later store initial values of the path matrix
        p_m = list()

        #iterates through the number of states to create an empty matrix of the
        #appropriate size for the score matrix
        for i in range(len(self.states)):
            s_m.append([0.0]*len(self.seq))
        #defining the score matrix as a data frame so it is more
        #intuitive to call values from the matrix
        score_matrix = pd.DataFrame(np.matrix(s_m), columns = pos, index = self.states)

        #iterate through the sequence and calculates the probability/score for
        #each state
        #score = sum of previous scores * emit(syb_i,state_i)
                    #* tran(state_i-1, state_i)
        for i in pos:
            #the first symbol in the sequence is treated differently because
            #it has not previous scores/probabilities to calculate its value
            if i == 0:
                for s in self.states:
                    #score for the first symbol =
                                #1 * 1/(#of states) * emit(syb_0,state_s)
                    #transition value = 1/(#of states)
                    score_matrix[i][s] = (1/len(self.states))*self.emission[self.seq[i]][s]
            else:
                #iterates through the states we could go into next
                for next in self.states:
                    #iterates through the states we could have come from
                    for prev in self.states:
                        #defining the values of each part of the
                        #probability/score calculation to make the code
                        #more readable
                        prev_score = score_matrix[i-1][prev]
                        tran = self.transition[next][prev]
                        emit = self.emission[self.seq[i]][next]
                        #calculating the probablity/score for the state
                        #and symbol in the matrix
                        score_matrix[i][next] += prev_score*tran*emit
        #calculate the probabilty of the emitted sequence by summing all of
        #the values in the matrix
        prob = float()
        for s in self.states:
            prob += score_matrix[len(self.seq)-1][s]
        return prob



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


    my_HMM = HMM(seq, states, transition, emission)
    #calculates the probability of the emitted sequence for our given HMM
    prob = my_HMM.calcProb()
    #prints the probability
    print(prob)



if __name__ == '__main__':
    main()
