#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-11-06
# Github acct: edevog
import sys
import numpy as np
import pandas as pd
'''
This program implements Baum-Welch learning in order to determine the
transition and emission probability matrices for an HMM given the emitted
sequence and the initial transition and emission probability matrices.
Input:
A text file with the number of iterations, the emitted sequence, the symbols
that compose the emitted sequence, the possible states, and the initial
emission and transition probability matrices
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
    def __init__(self, seq, states, transition, emission, sybs):
        self.seq = seq
        self.states = states
        self.sybs = sybs
        self.transition = transition
        self.emission = emission


    def forward(self):
        '''
        The forward algorithm
        '''
        #an empty list to later store initial values of the score matrix
        s_m = list()


        #iterates through the number of states to create an empty matrix of the
        #appropriate size for the score matrix
        for i in range(len(self.states)):
            s_m.append([0.0]*len(self.seq))
        #defining the score matrix as a data frame so it is more
        #intuitive to call values from the matrix
        forward_matrix = pd.DataFrame(np.matrix(s_m), columns = list(range(0, len(self.seq))), index = self.states)

        #iterate through the sequence and calculates the probability/score for
        #each state
        #score = sum of previous scores * emit(syb_i,state_i)
                    #* tran(state_i-1, state_i)
        for i in range(len(self.seq)):
            #the first symbol in the sequence is treated differently because
            #it has not previous scores/probabilities to calculate its value
            if i == 0:
                for s in self.states:
                    #score for the first symbol =
                                #1 * 1/(#of states) * emit(syb_0,state_s)
                    #transition value = 1/(#of states)
                    # print((1/len(self.states))*self.emission[self.seq[i]][s])
                    forward_matrix[i][s] = (1/len(self.states))*self.emission[self.seq[i]][s]
            else:
                #iterates through the states we could go into next
                for next in self.states:
                    #iterates through the states we could have come from
                    for prev in self.states:
                        #defining the values of each part of the
                        #probability/score calculation to make the code
                        #more readable
                        prev_score = forward_matrix[i-1][prev]
                        tran = self.transition[next][prev]
                        emit = self.emission[self.seq[i]][next]
                        #calculating the probablity/score for the state
                        #and symbol in the matrix
                        forward_matrix[i][next] += prev_score*tran*emit
        prob = float()
        for s in self.states:
            prob += forward_matrix[len(self.seq)-1][s]
        return forward_matrix, prob

    def backward(self):
        '''
        The backward algorithm
        '''
        b_tran = self.transition.T
        b_transition = pd.DataFrame(b_tran, columns = self.states, index = self.states)

        #an empty list to later store initial values of the score matrix
        s_m = list()
        #reverse the sequence
        b_seq = []
        for i in reversed(self.seq):
            b_seq.append(i)
        # iterates through the number of states to create an empty matrix of the
        # appropriate size for the score matrix
        for i in range(len(self.states)):
            s_m.append([0.0]*len(b_seq))
        #defining the score matrix as a data frame so it is more
        #intuitive to call values from the matrix
        backward_matrix = pd.DataFrame(np.matrix(s_m), columns = list(range(0, len(self.seq))), index = self.states)

        #iterate through the sequence and calculates the probability/score for
        #each state
        #score = sum of previous scores * emit(syb_i,state_i)
                    #* tran(state_i-1, state_i)
        for i in range(len(self.seq)):
            #the first symbol in the sequence is treated differently because
            #it has not previous scores/probabilities to calculate its value
            if i == 0:
                for s in self.states:
                    #score for the first symbol =
                                #1 * 1/(#of states) * emit(syb_0,state_s)
                    #transition value = 1/(#of states)
                    backward_matrix[i][s] = 1
            else:
                for next in self.states:
                    for prev in self.states:

                        prev_score = backward_matrix[i-1][prev]
                        transit = b_transition[next][prev]
                        emit = self.emission[b_seq[i-1]][prev]
                        backward_matrix[i][next] += prev_score*transit*emit
        prob = float()
        for s in self.states:
            prob += backward_matrix[len(self.seq)-1][s]
        return backward_matrix, prob

    def makeNodeResponsibilityMatrix(self, backward_matrix, forward_matrix, emission, f_sink):
        '''
        creates the node responsiblity matrix which determines the probability
        the HMM is in a particular state at a particular point in the emitted
        sequence for each state and emitted symbol
        '''
        prob_dict = {}
        for k in self.states:
            prob_dict[k] = []
        for i in range(len(self.seq)):
            for k in self.states:
                prob_dict[k].append(forward_matrix[i][k] * backward_matrix[len(self.seq)- 1 - i][k]/f_sink)
                prob = forward_matrix[i][k] * backward_matrix[len(self.seq)- 1 - i][k]/f_sink
        return prob_dict

    def makeEdgeResponsibilityMatrix(self, forward_matrix, transition, emission, backward_matrix, f_sink):
        '''
        creates the edge responsiblity matrix which determines the probability
        that the HMM passes through the edge connecting two states given the
        HMM emitted a specific symbol
        '''
        edge_matrix = {l : {k : [] for k in self.states} for l in self.states}
        for i in range(len(self.seq)-1):
            for l in self.states:
                forward_li = forward_matrix[i][l]
                for k in self.states:
                    backward_ki = backward_matrix[len(self.seq) - 2 - i][k]
                    weight_ilk = emission[self.seq[i+1]][k]*transition[k][l]
                    edge_matrix[l][k].append(forward_li*backward_ki*weight_ilk/f_sink)
        return edge_matrix


    def estimateEmissionMatrix(self, node_responsibility_matrix):
        '''
        determines the emission matrix from the node responsibility matrix
        by summing the probabilities then dividing by the total of each row
        '''
        e_m = np.zeros((len(self.states),len(self.sybs)),dtype = "longdouble")
        emission_matrix = pd.DataFrame(e_m, columns = self.sybs, index = self.states)
        norm_emit = pd.DataFrame(e_m, columns = self.sybs, index = self.states)
        for l in self.states:
            for s in self.sybs:
                for i in range(len(self.seq)):
                    if s == self.seq[i]:
                        emission_matrix[s][l] += node_responsibility_matrix[l][i]

        totals = []
        for l in self.states:
            total = float()
            for s in self.sybs:
                total += emission_matrix[s][l]
            totals.append(total)
        for i in range(len(self.states)):
            for s in self.sybs:
                norm_emit[s][self.states[i]] = emission_matrix[s][self.states[i]]/totals[i]
        self.emission = norm_emit
        return norm_emit

    def estimateTransitionMatrix(self, edge_responsibility_matrix):
        '''
        determines the transition matrix from the edge responsiblity matrix
        by summing the probabilities then dividing by the total of each row
        '''
        t_m = np.zeros((len(self.states),len(self.states)),dtype = "longdouble")
        transition_matrix = pd.DataFrame(t_m, columns = self.states, index = self.states)
        norm_transit = pd.DataFrame(t_m, columns = self.states, index = self.states)
        totals = []
        for l in self.states:
            for k in self.states:
                transition_matrix[l][k] = sum(edge_responsibility_matrix[k][l])
        for l in self.states:
            total = []
            for k in self.states:
                total.append(transition_matrix[k][l])
            totals.append(sum(total))
        for i in range(len(self.states)):
            for l in self.states:
                norm_transit[l][self.states[i]] = transition_matrix[l][self.states[i]]/totals[i]
        self.transition = norm_transit
        return norm_transit

    def calcEstimates(self, emission, transition):
        forward_matrix, f_sink = self.forward()
        backward_matrix, b_sink = self.backward()
        node_responsibility_matrix = self.makeNodeResponsibilityMatrix(backward_matrix, forward_matrix, emission, f_sink)
        edge_responsibility_matrix = self.makeEdgeResponsibilityMatrix(forward_matrix, transition, emission, backward_matrix, f_sink)
        emission = self.estimateEmissionMatrix(node_responsibility_matrix)
        transition = self.estimateTransitionMatrix(edge_responsibility_matrix)
        # emission = estimated_emission
        # transition = estimated_transition
        return emission, transition


def main():
    #this chunk of code parses the datafile
    iter = int(sys.stdin.readline().strip())
    #the rest of the data file contents
    f = sys.stdin.readlines()
    seq = list(f[1].strip())
    #the unique symbols that occur in the emitted sequence
    sybs = f[3].strip().split()
    #the possible states
    states = f[5].strip().split()
    #the matrix of transition probabilities
    tran_matrix = []
    #the for loop allows us to define matrices of different dimensions based
    #on the number of states in the HMM
    for i in range(len(states)):
        tran_matrix.append(list(map(float, f.pop(8).strip().split()[1:])))
    tran_matrix = np.array(tran_matrix, dtype = np.longdouble)
    #defines the transition matrix as a dataframe with the column and row
    #names as the states so it is more intuitive to call data from the matrix
    # b_tran = tran_matrix.T

    # b_transition = pd.DataFrame(b_tran, columns = states, index = states)

    transition = pd.DataFrame(tran_matrix, columns = states, index = states)
    #the matrix of emission probabilities
    emit_matrix = list()
    #the for loop allows us to define matrices of different dimensions based on
    #the number of states
    for i in range(len(states)):
        emit_matrix.append(list(map(float, f.pop(10).strip().split()[1:])))
    emit_matrix = np.array(emit_matrix, dtype = np.longdouble)
    #defines the emission matrix as a data frame with the column names as the
    #unique symbols and the row names as the states so it is more intuitive to
    #call data from the matrix
    emission = pd.DataFrame(emit_matrix, columns = sybs, index = states)

    my_HMM = HMM(seq, states, transition, emission, sybs)
    #runs the baum-welch algorithm
    for i in range(iter):
        emission, transition = my_HMM.calcEstimates(emission, transition)
        if i == iter-1:
            print(round(transition,3).to_csv(sep="\t"))
            print("--------")
            print(round(emission,3).to_csv(sep="\t"))
if __name__ == '__main__':
    main()
