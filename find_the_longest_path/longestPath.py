#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-10-31
# Github acct: edevog

import sys
import re

'''
This program finds the longest path through a directed graph given a source and
sink node. Each edge in the graph has a weight. Our longest path will traverse
edges so that the sum of all the edges' weights along our path is the largest.

Input: A text file where the frist two lines are the source node and sink node
respectively. The rest of the text file is an adjacency graph with edge weights
(node->out node:edge weight)
0
4
0->1:7
0->2:4
2->3:2
1->4:1
3->4:3

Output:
The first line is the path's score. The second line is the path taken delimited
by arrows.
9
0->2->3->4
'''

class LongestPath:
    '''
    Finds the longest path in the provided adjacency graph.

    Attributes:
        data = adjacency graph provided by the data file from standard in
        i_node = the source node specified in the data file
        o_node = the sink node specified in the data file
    '''
    def __init__(self, data, i_node, o_node):
        self.data = data
        self.i_node = i_node
        self.o_node = o_node

    def makeOutDict(self):
        '''
        creates a diction of nodes as keys, their out node, weight, and a boolean to
        track if they have been used as the values
        '''
        out_dict = {}
        for d in self.data:
            #for each line in the adjacency graph split it into a list
            #on the arrow and colon
            l = re.split('\->|:', d)
            #remove the first value of the list and store it as the node
            n = l[0]
            #remove the first value of the updated list and store it as the
            #out node
            out = l[1]
            #set the final value in the list as the weight and define it as
            #an integer
            wt = int(l[2])
            #if the node is not already in the dictionary, store it as a key
            #and set the output node and weight as its values
            if n not in out_dict:
                out_dict[n] = [[out], [wt]]
            else:
                #if the node is already in the dictionary, append the out node
                #and weight values to the exisiting list
                out_dict[n][0].append(out)
                out_dict[n][1].append(wt)
        return out_dict

    def makeInDict(self, out_dict):
        '''
        Makes an input dictionary using the existing dictionary
        of outputs (matchDict)
        '''
        in_dict = {}
        unique = []
        #iterate through each node in the output dictionary,
        for key in out_dict:
            #iterate through each output of the node
            for i in range(len(out_dict[key][0])):
                #set the input node as the ith output node
                n_in = out_dict[key][0][i]
                #set the input node's weight as the ith weight of the
                #output node
                wt = out_dict[key][1][i]
                #if the input node is not already in the dictionary,
                if n_in not in in_dict:
                    unique.append(n_in)
                    #create a new key as the input node and the values as the
                    #key from the output dictionary and the input node's weight
                    in_dict[n_in] = [[key], [wt]]
                #if the input node is already in the dictionary,
                else:
                    #append the key from the output dictionary and the
                    #input node's weight to the existing list
                    in_dict[n_in][0].append(key)
                    in_dict[n_in][1].append(wt)
            #after iterating through all the output nodes check if any nodes in
            #the output dictionary are not in the input dictionary and add
            #any missing nodes so the input dictionary contains all nodes in
            #the graph
            if key not in in_dict:
                in_dict[key] = [[],[]]
        return in_dict

    def makeIODicts(self):
        '''
        Calls the functions to make the input and output dictionaries, then
        updates the output dictionary to include all the nodes in the input
        dictionary so the output dictionary has every node in the graph
        '''
        #make the output dictionary
        out_dict = self.makeOutDict()
        #make the input dictionary
        in_dict = self.makeInDict(out_dict)
        #if any nodes in the input dictionary are not in the output dictionary,
        #add them to the output dictionary so it contains all nodes in the graph
        for key in in_dict:
            if key not in out_dict:
                out_dict[key] = [[],[]]
        return out_dict, in_dict

    def makeScoreDict(self, out_dict, in_dict):
        '''
        This dictionary stores the scores for each node as we find the
        longest path.
        '''
        score_dict = {}
        count_dict = self.makeCountDict(out_dict, in_dict)
        for n in count_dict:
            if n != self.i_node:
                #we want all source nodes to have a lower score than our
                #specified source node
                score_dict[n] = [-sys.maxsize, str()]
            elif n == self.i_node:
                #any nodes that have an input will have their score calculated
                #when we iterate through our topological order
                score_dict[n] = [int(), str()]

        return score_dict



    def makeCountDict(self, out_dict, in_dict):
        '''
        make a dictionary of nodes and a count of how many inputs they have
        '''
        count_dict = {}
        #iterates through each node in the output dictionary
        for n in out_dict:
            #if the node is not in the input dictionary, store it's count in
            #the count dictionary as 0
            if n not in in_dict:
                count_dict[n] = 0
            else:
                #if the node is in the output dictionary, store it's count
                #as the length of its input node list
                count_dict[n] = len(in_dict[n][0])
        #iterate through all nodes in the input dictionary
        for n in in_dict:
            # if the node is not in the count dictionary,
            if n not in count_dict:
                #store its count as the length of the input node list
                count_dict[n] = len(in_dict[n][0])
        return count_dict

    def makeTopOrd(self, count_dict, out_dict):
        '''
        Returns a list of the topological order of the graph determined using
        Kahn's algorithm
        '''
        q = []
        top_ord = []
        #iterate through the count dict to find all the nodes with 0 inputs
        for n in count_dict:
            if count_dict[n] == 0:
                #appending all of these nodes to our queue list (q)
                q.append(n)
        #while there are still items in the queue,
        while len(q) > 0:
            #we remove one of the nodes in the queue
            j = q.pop()
            #append it to our topological order list
            top_ord.append(j)
            #iterate through all the output nodes for our selected node
            for i in out_dict[j][0]:
                #subtract one from the count dictionary for the output node
                count_dict[i] += -1
                #check if the output node now has a count of 0
                if count_dict[i] == 0:
                    #add it to the queue
                    q.append(i)
        return top_ord


    def calcScore(self, node, in_dict, score_dict):
        '''
        calculates the score for a node and adds it to the score dictionary
        if it is better than the exisiting score
        '''
        #checking if the node has an input
        if len(in_dict[node][0]) != 0:
            #iterate through the length of the values in the in_dict for a node
            for i in range(len(in_dict[node][0])):
                #-------------------------------------------------------------
                #all of these variables are defined in order to make the logic
                #of the code more easily understood without having to look up
                #the structure of the dictionaries
                #-------------------------------------------------------------
                #the ith in_node for our node
                in_node = in_dict[node][0][i]
                #the ith weight of our in_node
                in_wt = in_dict[node][1][i]
                #the score of our ith in_node
                in_score = score_dict[in_node][0]
                #the ith in_node's score's path
                # in_path = score_dict[in_node][1]
                #the new score for our node
                new_score = in_score + in_wt
                #checking if our new score is better than the existing score in
                #the score dictionary
                if score_dict[node][0] < new_score:
                    new_path = in_node
                    score_dict[node] = [new_score, in_node]
        return score_dict

    def findLongestPath(self):
        '''
        Iterates through the topological order of our graph and calculates
        the score for each node. Then it returns the score and path of the
        longest path for our sink node
        '''
        out_dict, in_dict = self.makeIODicts()
        count_dict = self.makeCountDict(out_dict, in_dict)
        top_ord = self.makeTopOrd(count_dict, out_dict)
        score_dict = self.makeScoreDict(out_dict, in_dict)
        node = self.o_node
        path= [self.o_node]
        i=0
        for n in top_ord:
            i+=1
            score_dict = self.calcScore(n, in_dict, score_dict)
        #defines the score of the longest path for our sink node
        score = score_dict[self.o_node][0]
        #defines the path of the longest path

        while len(score_dict[node][1]) != 0:
            path.append(score_dict[node][1])
            node = score_dict[node][1]
        path.reverse()

        return score, path





def main():
    i_node = sys.stdin.readline().strip()
    o_node = sys.stdin.readline().strip()
    data = [d.strip() for d in sys.stdin.readlines()]

    longest_path = LongestPath(data, i_node, o_node)

    score, path = longest_path.findLongestPath()

    print(score)
    print("->".join(path))

if __name__ == '__main__':
    main()
