#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-10-27
# Github acct: edevog

import sys
import random
import re

'''
This program constructs a string from kmers.

Input: Text file with kmer length and all kmers
4
CTTA
ACCA
TACC
GGCT
GCTT
TTAC

Output:
A string of all connected kmers (GGCTTACCA). There are potentially many strings
that could be constructed. This program only outputs one of the potential
strings.
'''

class DeBruijn(object):
    '''
    The functions necessary to create a De Bruijn graph.

    Attributes:
        kmers = the list of kmers from the input data file.
        matchDict = a dictionary of each kmer node and its outputs
    '''
    def __init__(self, kmers):
        self.ks = kmers
        self.matchDict = {}

    def findMatch(self):
        '''
        Creates nodes for each kmer in the provided data file. Then finds the
        nodes that "match" each other --- the node whose first letters
        matches the last letters of another node.
        '''
        nodes = []
        for kmer in self.ks:
            #creates a list of nodes
            nodes.append(kmer[:-1])
            nodes.append(kmer[1:])
        #sorts the nodes in alphabetical order
        srt_nodes = sorted(nodes)
        for n in srt_nodes:
            #creates a dictionary with nodes as keys and values of an
            #empty list
            self.matchDict.update({n:[]})
        for k in self.ks:
            #finds the matching nodes for each node key and appends it to its
            #value list
            for key in self.matchDict.keys():
                if k[:-1] == key:
                    self.matchDict[key].append(k[1:])

class Node(object):
    '''
    All functions required to modify and find the attributes of a node

    Attributes:
        matchDict = a dictionary of every kmer node and its outputs
        inputDict = a dictionary of every kmer node and its inputs
        nodes = a list of nodes created from the provided kmers
        usedOutputs = a list of outputs from a node that have already been
                      traversed
        node_o = the initial node our string will start with
        node_n = the terminal node our string will end with
    '''
    def __init__(self, matchDict, inputDict):
        self.matchDict = matchDict
        self.inputDict = inputDict
        self.nodes = []
        self.usedOutputs = []
        self.node_o = ""
        self.node_n = ""

    def findInputs(self, name):
        '''
        Finds the inputs of a node by referencing the inputDict
        '''
        inputs = self.inputDict[name]
        return inputs

    def findOutputs(self, name):
        '''
        Finds the outputs of a node by referencing the matchDict
        '''
        outputs = self.matchDict[name]
        return outputs

    def freeOutputs(self, name, usedNodeIndex):
        '''
        Removes used outputs of a node
        '''
        outputs = self.matchDict[name]
        outputs.pop(usedNodeIndex)
        return outputs


class FindSeq:
    '''
    Creates an Eularian cycle in a graph from a De Bruijn graph. Finds a path
    through the Eularian cycle that uses every output from every node. Joins
    each node in the path into a sequence.

    Attributes:
        matchDict = a dictionary of all kmer nodes and their outputs
        inputDict = an empty dictionary that is populated in this class with
                    every node and its inputs

    '''

    def __init__(self, matchDict):
        self.matchDict = matchDict
        self.inputDict = {}


    def makeDicts(self):
        '''
        creates the output dictionary (matchDict)
        calls the function to make the input dictionary (inputDict)
        connects the initial node and terminal node to create an Eularian cycle
        '''
        self.inputDict = self.makeInputDict()
        #checks that the inputDict and matchDict have the same keys
        for key in self.inputDict:
            if key not in self.matchDict:
                self.matchDict[key] = []
        for key in self.matchDict:
            if key not in self.inputDict:
                self.inputDict[key] = []
        #finds the initial and terminal nodes
        self.node_o, self.node_n = self.findNodeEnds()
        node_o, node_n = self.findNodeEnds()
        #adds the terminal node as an input node for the initial node
        if node_o in self.inputDict:
            self.inputDict[node_o].append(node_n)
        else:
            self.inputDict[node_o] = [node_n]
        #adds the initial node as an output for the terminal node
        if node_n in self.matchDict:
            self.matchDict[node_n].append(node_o)
        else:
            self.matchDict[node_n] = [node_o]
        #now we have two dictionaries containing all the info for an
        #Eularian *Cycle*
        return self.matchDict, self.inputDict

    def makeInputDict(self):
        '''
        Makes an input dictionary using the existing dictionary
        of outputs (matchDict)
        '''
        self.inputDict = {}
        for key in self.matchDict:
            #creates a key for each value (output node) in the matchDict
            for value in self.matchDict[key]:
                if value not in self.inputDict:
                    self.inputDict[value] = []
                #adds the key from the matchDict as a value in the inputDict
                self.inputDict[value].append(key)
        return self.inputDict

    def findNodeEnds(self):
        '''
        Identifies the inital node and the terminal node
        '''
        node_o = ""
        node_n = ""

        for key in self.matchDict.keys():
            if len(self.matchDict[key]) < len(self.inputDict[key]):
                node_n = key
            if len(self.matchDict[key]) > len(self.inputDict[key]):
                node_o = key
        return node_o, node_n


    def randomNode(self):
        '''randomly selects a node'''
        random.choice(self.matchDict.keys())

    def traverseNodes(self, i_node, path):
        '''
        Goes through each node's outputs to construct the string of all
        connected nodes.
        '''
        #takes the incoming list and makes a new list so it can append new
        #values to the incoming list
        #without this step, the values are appended in a list of lists
        #(ie path = [[n1, n2, n3], n4, n5])
        flat_path = []
        for i in path:
            if len(i) > 1:
                for j in i:
                    flat_path.append(j)
        else:
            flat_path = path
        myNode = Node(self.matchDict, self.inputDict)
        curr_node = i_node
        outputs = myNode.findOutputs(curr_node)
        #checks if there are any outputs that have not been traversed for
        #a node
        while len(outputs) > 0:
            #randomly selects an output from the list of available outputs
            #and sets it as n_node
            x = random.choice(range(len(myNode.findOutputs(curr_node))))
            n_node = myNode.findOutputs(curr_node)[x]
            #appends the node to our existing path list
            flat_path.append(n_node)
            outputs = myNode.freeOutputs(curr_node, x)
            n_outputs = myNode.findOutputs(n_node)
            outputs = n_outputs
            curr_node = n_node
        return flat_path


    def newPath(self, curr_path):
        '''
        Creates a new path that starts with the node with unused outputs
        '''
        myNode = Node(self.matchDict)
        for node in curr_path:
            if len(self.matchDict[node]) > 0:
                i = curr_path.index(node)

                n_path = curr_path[i:] + curr_path[1:i+1]
                return n_path


    def findPath(self, r_node, curr_path):
        '''
        Iterates through the nodes until all output nodes are included in
        the path
        '''
        values_list = []
        count = 0
        n_path = [r_node]
        i_node = r_node
        v_list_len = sum(map(len, self.matchDict.values()))
        while len(curr_path) < v_list_len:
            curr_path = self.traverseNodes(i_node, n_path)
            #if the current path does not contain all edges
            #find a node on the path with unused outputs and create a new
            #path starting at that node
            x_path = []
            if len(curr_path) < v_list_len:
                for node in curr_path:
                    if len(self.matchDict[node]) > 0:
                        i = curr_path.index(node)
                        #create a new path that starts with the node with
                        #unused outputs
                        x_path.append(curr_path[i:] + curr_path[1:i+1])
                n_path = x_path[0]
                i_node = n_path.pop()
                n_path.append(i_node)
            else:
                return curr_path
    def outPath(self, f_path):
        '''
        Takes the final path and splits it so it starts with the initial node
        and ends with the terminal node. Joins each node together to form a
        string sequence
        '''
        n_index = int()
        match = []
        oddList = []

        for i in range(len(f_path)-1):
            if f_path[i] == self.node_n and f_path[i+1] == self.node_o:
                n_index = i
        # splits the final path so it starts with the initial node and ends
        #with the terminal node
        i_path = f_path[n_index+1:]+f_path[1:n_index+1]

        #joins each node to the previous node
        for kmer in i_path:
            for j in range(len(i_path)):
                if kmer[1:] == i_path[j][:-1]:
                    match.append(kmer)
                    match.append(i_path[j])

        #removes duplicate nodes, which occur every other time
        a = match[0]
        allList = list(range(len(match)))
        for i in allList:
            if i%2 != 0:
                oddList.append(i)
        #starts the sequence with the fist nodes from the match list
        seq = [a]
        len_kmer = len(a)
        #appends each node with the last character from the next node since the
        #middle characters overlap
        for i in oddList:
            seq.append(match[i][len_kmer-1:])
        return("".join(seq))



def main():
    #defines the first line of the file as len_kmer
    len_kmer = int(sys.stdin.readline())
    #defines the other lines of the file as kmers
    data = sys.stdin.readlines()
    kmers = [d.strip() for d in data]

    #makes a De Bruijn graph from the kmers
    myDeBruijn = DeBruijn(kmers)
    matches = myDeBruijn.findMatch()


    #finds a path which traverses every node, starts with the initial node and
    #ends on the terminal node
    myfindSeq = FindSeq(myDeBruijn.matchDict)
    #populates the dictionaries within the FindSeq class
    outputDict, inputDict = myfindSeq.makeDicts()
    r_node = random.choice(list(myfindSeq.matchDict.keys()))
    path = [r_node]
    f_path = myfindSeq.findPath(r_node, path)
    seq = myfindSeq.outPath(f_path)
    #prints the sequence comprized of our given kmers as a string
    print(seq)



if __name__ == '__main__':
    main()
