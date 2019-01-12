#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-12-04
# Github acct: edevog
import sys

'''
This program determines the cyclospectrum of a given peptide.
Input: A text file of a peptide
Ex:
LEQN

Output:
The cyclospectrum of the peptide in a space delimited list.
Ex:
0 113 114 128 129 227 242 242 257 355 356 370 371 484
'''

class Peptide:
    """
    This class contains methods to determine the properties of a peptide.
    Attributes:
        peptide = the given peptide
        aa_mass = a dictionary of the amino acids and their mass
    """
    def __init__(self, peptide):
        self.peptide = peptide
        self.aa_mass = self.makeAminoAcidMassDict()

    def makeAminoAcidMassDict(self):
        '''
        a dictionary of the interger masses of each amino acid
        '''
        aa_mass = {
        "G" : 57,
        "A" : 71,
        "S" : 87,
        "P" : 97,
        "V" : 99,
        "T" : 101,
        "C" : 103,
        "I" : 113,
        "L" : 113,
        "N" : 114,
        "D" : 115,
        "K" : 128,
        "Q" : 128,
        "E" : 129,
        "M" : 131,
        "H" : 137,
        "F" : 147,
        "R" : 156,
        "Y" : 163,
        "W" : 186
        }
        return aa_mass

    def makeSubPeptides(self):
        '''
        creates all subpeptides for a cyclic peptide
        '''
        c_peptide = self.peptide + self.peptide
        sub_peptides = []
        for i in range(len(self.peptide)):
            for j in range(len(self.peptide)):
                if i < i + j:
                    sub_peptides.append(c_peptide[i:j+i])
        sub_peptides.append(self.peptide)
        return sub_peptides

    def findMass(self, sub_peptides):
        '''
        finds the mass of each peptide in the subpeptide list
        '''
        spectrum = [0]
        for p in sub_peptides:
            total = 0
            for l in p:
                total += self.aa_mass[l]
            spectrum.append(total)
        return sorted(spectrum)


def main():
    #the first list of data file is the given peptide
    peptide = list(sys.stdin.readline().strip())
    my_peptide = Peptide(peptide)
    sub_peptides = my_peptide.makeSubPeptides()
    spectrum = my_peptide.findMass(sub_peptides)
    print(*spectrum, sep = ' ')

if __name__ == '__main__':
    main()
