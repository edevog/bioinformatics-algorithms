#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-10-10
# Github acct: edevog
import random
from fastaReader import *
import math

'''
This program finds the promoter motif from a list of sequences. We assume that
the promoter motif is the motif from the motif matrix with the lowest entropy.

Input:
A fasta file that presents the upstream 50 bases of a speicies
(ie Pyrobaculum group). For the fastareader to work, the first line of the
fasta file must have a carrot (>).

The command line argument should follow the format:
randomizedMotifSearch.py -i=100000 -p=1 -k=13 <somefile.fa >someOutputFIle.fa

where i is the number of iterations, p is the number of psuedocounts, and
k is the motif length.

Output:
The promoter motif (the motif from the motif matrix with the lowest entropy)
and the motif entropy score

'''

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.

        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Program prolog - reads in a fasta file',
                                             epilog = 'Program epilog - file must begin with a carrot (>)',
                                             add_help = True,
                                             prefix_chars = '-',
                                             usage = '%(prog)s [p, i, k] -pseudoCounts[null] -iterations[null] -motifLength[null] <fastafile >txtfile'
                                             )
        self.parser.add_argument('-p', '--pseudoCounts', type = int, action = 'store', help='pseudocounts')
        self.parser.add_argument('-i', '--iterations', type = int, action = 'store', help='number of iterations')
        self.parser.add_argument('-k', '--motifLength', type = int, action = 'store', help='the length of the motif')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and
    eventual exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg



class FindMotif(object):
    """This class find the promoter motif."""
    def __init__(self, iterations, kmerLength, pseudoCounts):
        self.iterations = iterations
        self.kmerLength = kmerLength
        self.pseudoCounts = pseudoCounts
        #I am definiting sequences as a member variable because we access the
        #sequences from our dataset several times throughout this class
        self.sequences = []

        reader = FastAreader()
        for head, seq in reader.readFasta():
            self.sequences.append(seq)
        #I am defining kmerMetrix as a member variable because the matrix is



    def randomMotif(self, seq):
        '''
        picks  random kmer from each sequence
        '''
        seqLength = len(seq)
        listOfChoices = range(0,(seqLength-self.kmerLength))
        start = random.choice(listOfChoices)
        end = start + self.kmerLength
        randomMotif = seq[start:end]
        return randomMotif


    def createCountMatrix(self, motifMatrix):
        '''
        inputs pseudocounts then counts how many times each nucleotide base
        occurs at each position in the kmers. then divides by the kmer length
        '''
        #defining the bases found within our sequences
        bases = ['A', 'C', 'T', 'G']
        #defining a list of pseudocounts that are independent so when we add
        #to one list, we do not add to all of the lists
        pseudoCountList = [[self.pseudoCounts]*self.kmerLength, [self.pseudoCounts]*self.kmerLength, [self.pseudoCounts]*self.kmerLength, [self.pseudoCounts]*self.kmerLength]
        #creating a count matrix as a dictionary where the bases are the keys
        #and the psuedocounts are the values
        countMatrix = dict((zip(bases, pseudoCountList)))
        #iterate through the lenght of our defined kmer length
        for i in range(self.kmerLength):
            #iterate through the kmers in our motif matrix
            for kmer in motifMatrix:
                #each time we encounter the ith base in our kmer, add one to
                #the ith base in our count matrix in the ith column of the
                #count matrix
                countMatrix[kmer[i]][i] += 1
                #in other words we are iterating through the base at every
                #position of each kmer in our matrix then tracking how many
                #times that base occurs in that position and storing that
                #information in the count dict
        #iterate through the list of counts for each base
        for countList in countMatrix.values():
            #iterate through each position's count in the list, divde the
            #count by the number of sequences then add the psuedocounts
            #value multiplied by the number of bases (4)
            countList[:] = [x/(len(self.sequences)+(self.pseudoCounts*4)) for x in countList]
        return countMatrix

    def createConsensusMotif(self, countMatrix):
        '''
        creates the consensus motif base on the alues from the count matrix
        '''
        conMotif = []
        for i in range(0, self.kmerLength):
            score = 0
            bestNuc = ""
            for key in countMatrix:
                if countMatrix[key][i] > score:
                    score = countMatrix[key][i]
                    bestNuc = key
            conMotif.append(bestNuc)
        return "".join(conMotif)


    def scoreCountMatrix(self, countMatrix):
        '''
        calculates the entropy of the motif matrix by using the values from
        the countMatrix
        '''
        allColEntropies = []
        for i in range(self.kmerLength):
            colEntropies = []
            for base in countMatrix.keys():
                p = countMatrix[base][i]

                e = p*math.log(p, 2)
                colEntropies.append(e)
            colEntropy = sum(colEntropies)
            allColEntropies.append(colEntropy)
        entropy = sum(allColEntropies)*-1
        return entropy

    def scoreKmer(self, countMatrix):
        '''
        calculates the score using the probabilities int he countMatrix for
        each kmer
        '''
        listOfBestKmers = []
        for sequence in self.sequences:
            scores = []
            kmers = []
            #iterate through the number of kmers
            for i in range(len(sequence) - self.kmerLength+1):
                #defines the kmer we are working with
                kmer = sequence[i:i+self.kmerLength]
                kmers.append(kmer)
                allProbabilityOfBase = []
                #iterates through the position of bases in a kmer
                for j in range(self.kmerLength):
                    #defines our base
                    kmerBase = kmer[j]
                    #finds the probability value of the base
                    probabilityOfBase = countMatrix[kmerBase][j]
                    #stores it in a list
                    allProbabilityOfBase.append(probabilityOfBase)
                #initialize the score as 1 because we are multiplying
                score = 1
                #iterate through our list of probabilities
                for x in allProbabilityOfBase:
                    #multiply all the probabilities together to get a score
                    score *= x
                #store the score in a list
                scores.append(score)
            #finding the index of the best score in the list
            indexOfBestKmer = scores.index(max(scores))
            #using the index to pull the best kmer in the list
            bestKmer = kmers[indexOfBestKmer]
            #store all best kmers in a list
            listOfBestKmers.append(bestKmer)
        return listOfBestKmers


    def bestKmers(self):
        '''
        finds the best kmer from our initial random starting point
        '''

        count = 0
        randomMotifList = []
        for seq in self.sequences:
            slices = self.randomMotif(seq)
            randomMotifList.append(slices)
        motif_o = randomMotifList


        countMatrix_o = self.createCountMatrix(motif_o)
        #creates a matrix with the nucleotide bases and their counts for each position in a kmer
        scoreCountMatrix = self.scoreCountMatrix(countMatrix_o)
        consensusMotif_o = self.createConsensusMotif(countMatrix_o)

        #start a whie loop to continue until the initial score is less than the new score
        while True:
            #calculate the entropy of our initial count matrix
            consensusMotifScore_o = self.scoreCountMatrix(countMatrix_o)
            #pass the consensus motif over the sequences and create a new kmer matrix
            motif_new = self.scoreKmer(countMatrix_o)
            #calculate the counts for the new kmer matrix
            countMatrix_new = self.createCountMatrix(motif_new)
            #calculate the entropy of our new count matrix
            consensusMotifScore_new = self.scoreCountMatrix(countMatrix_new)
            #create a new consensus motif based on our new counts
            consensusMotif_new = self.createConsensusMotif(countMatrix_new)
            if consensusMotifScore_new >= consensusMotifScore_o:
                return consensusMotif_o, consensusMotifScore_o
            else:
                consensusMotif_o = consensusMotif_new
                countMatrix_o = countMatrix_new



def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere
    in process.

    '''
    if myCommandLine is None:
        myCommandLine = CommandLine()  # read options from the command line
    else :
        myCommandLine = CommandLine(myCommandLine) # interpret the list passed from the caller of main as the commandline.

    try:
        iterations = myCommandLine.args.iterations
        pseudoCounts = myCommandLine.args.pseudoCounts
        kmerLength = myCommandLine.args.motifLength

        ourMotif = FindMotif(iterations, kmerLength, pseudoCounts)

        bestKmers = []
        for i in range(iterations):
            bestKmer = ourMotif.bestKmers()
            bestKmers.append(bestKmer)

        bestKmerDict = {}
        for j in range(len(bestKmers)):
            bestKmerDict[bestKmers[j][0]] = bestKmers[j][1]

        winner = math.inf
        winnerKmer = ''
        for k in bestKmerDict.keys():
            if bestKmerDict[k] < winner:
                winner = bestKmerDict[k]
                winnerKmer = k

        print(winnerKmer, bestKmerDict[winnerKmer])



    except Usage as err:
       print (err.msg)

if __name__ == "__main__":
    main()
