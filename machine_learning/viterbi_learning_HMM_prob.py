#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-11-24
# Github acct: edevog
import sys
import pandas as pd
import itertools
import numpy as np

'''
This program determines the emission and transition probability matricies of
an HMM given an emitted sequence and the initial emission and transition
probability matricies using the Viterbi Learning algorithm.
Input:
A text file with the number of iterations to complete, the emitted sequence,
the symbols that comprise the emitted sequence, the states the emitted
sequence can be in, and the transition and emission matricies.
Ex.
100
--------
xxxzyzzxxzxyzxzxyxxzyzyzyyyyzzxxxzzxzyzzzxyxzzzxyzzxxxxzzzxyyxzzzzzyzzzxxzzxxxyxyzzyxzxxxyxzyxxyzyxz
--------
x   y   z
--------
A   B
--------
    A   B
A   0.582   0.418
B   0.272   0.728
--------
    x   y   z
A   0.129   0.35    0.52
B   0.422   0.151   0.426

Output:
A transition and emission probability matrix of the provided emitted sequence.
Ex.
A   B
A   0.875   0.125
B   0.011   0.989
--------
    x   y   z
A   0.0 0.75    0.25
B   0.402   0.174   0.424
'''

class HMM:
    '''
    Calculates and defines values related to a Hidden Markov Model

    Attributes:
        seq = the emitted sequence from the HMM
        states = the states in the HMM
        transition = the dataframe of the probabilities of transitioning
                        between states
        emission = the dataframe of the probabilities of being in a certain
                    state for a certain symbol in the sequence
    '''
    def __init__(self, seq, states):
        self.seq = seq
        self.states = states

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
                    path_matrix[i][s] = s
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
                            path_matrix[i][next] = path_matrix[i-1][prev] + next
        f_score = 0
        f_path = ""
        #iterate through the states in the final sequence value
        for s in states:
            #store the larger score and its path
            if f_score < score_matrix[len(seq)-1][s]:
                f_score = score_matrix[len(seq)-1][s]
                f_path = path_matrix[len(seq)-1][s]
        return f_path

    def tran_prob(self, states, path):
        '''
        This method determines the transition probability matrix based on the
        hidden path provided in the data file. First, it creates every
        possible pair that could occur from the given states. Next, it counts
        how many times each state transitions to another state; in other words
        it counts the number of Ax's, Bx's, and Cx's. If one of the states
        never occurs in the path, it sets the count equal to the number of
        states and sets the number of transition count to 1. Then, it
        counts the number of times each pair occurs in the path. Finally, it
        calculates the probability of each transition by dividing the pair
        counts by the transition count and stores it in an array labeled with
        the states.
        '''
        #initialize an empty dictionary to store the pair counts and the
        #state transition counts
        # tran_dict = {}
        tran_matrix = []
        # iterates through a list of tuples of all combinations
        tran_dict = {"".join(k):[0,0] for k in itertools.product(states, repeat = 2)}

        for i in range(len(path)-1):
            if path[i] not in tran_dict:
                tran_dict[path[i]] = 1
            else:
                tran_dict[path[i]] += 1
            tran_dict[path[i]+path[i+1]][0] += 1
        for p in tran_dict:
            if len(p) > 1:
                if p[0] not in tran_dict:
                    tran_dict[p] = [1,len(states)]
                else:
                    tran_dict[p][1] = tran_dict[p[0]]

        #determines the probability by diving the pair counts by the
        #transition count
        for s in states:
            row = []
            for p in states:
                if tran_dict[s+p][0] == 0:
                    row.append(0)
                else:
                    row.append(tran_dict[s+p][0]/tran_dict[s+p][1])
            tran_matrix.append(row)
        #stores the data in a dataframe to access intuitively later
        return pd.DataFrame(np.matrix(tran_matrix), columns = states, index = states)

    def emit_prob(self, seq, path, states, sybs):
        '''
        calculates the estimated emission probability
        '''
        emit_matrix = []
        emit_dict = {"".join(k):[0,0] for k in itertools.product(states, sybs)}

        for i in range(len(self.seq)):
            if path[i] not in emit_dict:
                emit_dict[path[i]] = 1
            else:
                emit_dict[path[i]] += 1
            emit_dict[path[i]+seq[i]][0] += 1
        for p in emit_dict:
            if len(p) > 1:
                if p[0] not in emit_dict:
                    emit_dict[p] = [1,len(states)]
                else:
                    emit_dict[p][1] = emit_dict[p[0]]

        #uses the above counts to calculate the emission probabilities
        for s in states:
            row = []
            for b in sybs:
                #if the symbol and state pair never happens the probability
                #of emission is 0
                if emit_dict[s+b][0] == 0:
                    row.append(0)
                #calculates the emission probability by dividing the count of
                #how many times the symbol and state occur in our given sequence
                #and path by the count of how many times the hidden path is in
                #the state
                else:
                    row.append(emit_dict[s+b][0]/emit_dict[s+b][1])
            emit_matrix.append(row)
        #stores the data in a dataframe to easily access later
        return pd.DataFrame(np.matrix(emit_matrix), columns = sybs, index = states)

    def findParams(self, sybs, iter, transition, emission):
        '''
        determines the hidden path from the emission and transition probability
        matrices, then determines the emission and transition probability
        matrices from the hidden path for the number of iterations specified
        in the data file
        '''
        i = 0

        while i < iter:
            i += 1
            path = self.findPath(self.seq, self.states, transition, emission)
            emission = self.emit_prob(self.seq, path, self.states, sybs)
            transition = self.tran_prob(self.states, path)
        return transition, emission


def main():
    #this chunk of code parses the datafile
    #the number of iterations we will use
    iter = int(sys.stdin.readline().strip())
    #the rest of the data file contents
    f = sys.stdin.readlines()
    #defines the emitted sequence from the HMM
    seq = list(f[1].strip())
    #the unique symbols that occur in the emitted sequence
    sybs = f[3].strip().split("\t")
    #the possible states
    states = f[5].strip().split("\t")
    #the matrix of transition probabilities
    tran_matrix = list()
    #the for loop allows us to define matrices of different dimensions based
    #on the number of states in the HMM
    for i in range(len(states)):
        tran_matrix.append(list(map(float, f.pop(8).strip().split()[1:])))
    #defines the transition matrix as a dataframe with the column and row
    #names as the states so it is more intuitive to call data from the matrix
    transition = pd.DataFrame(np.array(tran_matrix), columns = states, index = states)
    #the matrix of emission probabilities
    emit_matrix = list()
    #the for loop allows us to define matrices of different dimensions based on
    #the number of states
    for i in range(len(states)):
        emit_matrix.append(list(map(float, f.pop(10).strip().split()[1:])))
    #defines the emission matrix as a data frame with the column names as the
    #unique symbols and the row names as the states so it is more intuitive to
    #call data from the matrix
    emission = pd.DataFrame(np.array(emit_matrix), columns = sybs, index = states)

    my_HMM = HMM(seq, states)
    #calculates the new transition and emission probability matrices
    transition, emission = my_HMM.findParams(sybs, iter, transition, emission)
    print(round(transition,3).to_csv(sep = "\t"))
    print("--------")
    print(round(emission,3).to_csv(sep = "\t"))


if __name__ == '__main__':
    main()
