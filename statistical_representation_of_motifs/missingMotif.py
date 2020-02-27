
#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-10-02
# Github acct: edevog
import itertools as it
from fastaReader import *

'''
This program reads a fasta file from STDIN and ranks motifs based on how
statistically underrepresented the specific motif is.

Input:
A fasta file of a sequence where the first line contains a carrot (>) so the
fastareader works.

To call the program in the command line follow the following format:
    missingMotif.py -m 3 -M 8 -c -5

where m is the min motif length, M is the max motif length, and c is the
z-score cut off value

Output:
A report with headers ranking motifs by their kmer length
---longest to shortest--- then zscore---lowest to highest.
Includes the sequence with its reverse complement, its count, expected count,
and z-score

Ex:
sequence: reverse     count      Expect Zscore
AAAATTTA:TAAATTTT	835	1737.66	-21.66
GATTAATA:TATTAATC	550	1326.89	-21.33
AATTAATA:TATTAATT	929	1839.72	-21.24
AATTAATC:GATTAATT	977	1861.79	-20.51
.
.
.

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
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s [m, M, c] -minMotif[3] -maxMotif[8] -cutoff[-5] <fastafile >txtfile'
                                             )
        self.parser.add_argument('-m', '--minMotif', type = int, action = 'store', default = 3, help = 'this is the min kmer size')
        self.parser.add_argument('-M', '--maxMotif', type = int, action = 'store', default = 8, help = 'this is the max kmer size')
        self.parser.add_argument('-c', '--cutoff', type = float, action = 'store', default = -5, help = 'this is the Z-Score cut off')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)


class Usage(Exception):
    '''
    Used to signal a Usage error, evoking a usage statement and eventual
    exit when raised.
    '''
    def __init__(self, msg):
        self.msg = msg




class KmerData(object):
    """
    This class determines the expected count and z-score based on the actual
    count for every possible kmer

    Attributes:
        kmerMin = the minimum motif size to evaluate
        kmerMax = the maximum motif size to evaluate
        zScoreCutOff = the z-score cutoff
        kmerDict = a dictionary of all possible kmers for the given kmerMin and
                    kmerMax lengths provided and bases A, C, T, G as the keys
                    and the counts of the kmers as the values
    """
    def __init__(self, kmerMin, kmerMax, zScoreCutOff):
        self.kmerMin = kmerMin
        self.kmerMax = kmerMax
        self.zScoreCutOff = zScoreCutOff
        self.kmerDict = {}


        #creating every possible kmer using itertools and populated the keys
        #of a dictionary with the kmers
        for k in range(self.kmerMin, self.kmerMax + 1):
            for i in it.product('ACGT', repeat = k):
                self.kmerDict[''.join(i)] = [0]
            for kmer in self.kmerDict:
                if kmer == self.reverseComplement(kmer):
                    pass
                if kmer > self.reverseComplement(kmer):
                    self.kmerDict[kmer] = self.kmerDict[self.reverseComplement(kmer)]


    def reverseComplement(self, dna):
        #stackoverflow
        '''
        determines the reverse kmer of a sequence
        '''
        complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        try:
            return ''.join([complement[base] for base in dna[::-1]])
        except KeyError:
            pass

    def countKmers(self, seq):
        '''
        reads in all the kmers from the sequence, then counts how many times the kmer
        from the sequence matches a kmer in the dictionary
        '''
        for k in range(self.kmerMin, self.kmerMax + 1):
            for i in range(len(seq) - k + 1):
                kmerFromSeq = (seq[i:i+k])
                try:
                    if kmerFromSeq in self.kmerDict.keys():
                        self.kmerDict[kmerFromSeq][0] += 1
                except KeyError:
                    pass

    def calcExpectedCount(self, kmer):
        '''
        calculates the expected count of a kmer
        '''
        expectedCount = 0
        prek = self.kmerDict[kmer[:-1]][0]
        sufk = self.kmerDict[kmer[1:]][0]
        midk = self.kmerDict[kmer[1:-1]][0]
        expectedCount = (prek*sufk)/(midk)
        return expectedCount


    def calcZscore(self, kmer):
        '''
        calculates the Z-Score:
            zScore = (actual count - expected count)/sqrt(expected count)
        '''

        expectedCount = self.calcExpectedCount(kmer)
        kmerCount = self.kmerDict[kmer][0]
        zScore = 0
        zScore = (kmerCount - expectedCount)/(expectedCount**(0.5))
        return zScore

    def getData(self):
        '''
        calls the functions to calulate the expected count and zscore. appends
        the kmer, reverse kmer, actual count, expected count, and zscore values
        to a list of list
        '''
        outArray = []
        for kmer, count in self.kmerDict.items():
            if len(kmer) >= self.kmerMin + 2:
                kmerRC = self.reverseComplement(kmer)
                kmerWithComplements = kmer + ":" + kmerRC
                countVal = count[0]
                expCount = self.calcExpectedCount(kmer)
                zScore = self.calcZscore(kmer)
                if zScore <=  self.zScoreCutOff and kmer <= kmerRC:
                    outArray.append([kmer, kmerRC, count[0], self.calcExpectedCount(kmer), self.calcZscore(kmer)])
        return outArray


def main(myCommandLine=None):
    '''
    Implement the Usage exception handler that can be raised from anywhere in process.
    Defines the variables in the command line
    calls the classes
    prints the header
    prints the kmer, reverse kmer, count, expected count, and zscore

    '''
    if myCommandLine is None:
        # read options from the command line
        myCommandLine = CommandLine()
    else :
        # interpret the list passed from the caller of main as the commandline.
        myCommandLine = CommandLine(myCommandLine)

    try:
        minKmer = myCommandLine.args.minMotif - 2
        maxKmer = myCommandLine.args.maxMotif
        zScoreCutOff = myCommandLine.args.cutoff
        reader = FastAreader()
        kmerData = KmerData(minKmer, maxKmer, zScoreCutOff)

        for head, seq in reader.readFasta():
            kmerData.countKmers(seq)
        print("{0:8}:{1:8}\t{2:4}\t{3:7}\t{4:5}".format("sequence", "reverse", "count", 'Expect', "Zscore"))
        allData = kmerData.getData()
        #sorts data by length of sequence then by zscore
        sortedData = sorted(sorted(allData, key=lambda x: x[4]), key=lambda y: len(y[0]), reverse = True)

        for x in sortedData:
            print("{0:8}:{1:8}\t{2:0d}\t{3:.2f}\t{4:0.2f}".format(x[0], x[1], x[2], x[3], x[4]))

    except Usage as err:
       print (err.msg)


if __name__ == "__main__":
    main()
