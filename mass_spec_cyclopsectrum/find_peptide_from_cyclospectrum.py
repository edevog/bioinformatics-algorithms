#!/usr/bin/env python3
# Name: Beth DeVogelaere
# Date: 2018-12-04
# Github acct: edevog
import sys

'''
This program finds a cyclic peptide whose theoretical spectrum matches the
experimental spectrum given an ideal experimental spectrum.

Input:
A text file of a given ideal spectrum in a space delimited list.
Ex:
0 113 128 186 241 299 314 427

Output:
The unique mass profiles of the peptides whose cyclospectrum is within the
spectrum.
Ex:
186-128-113 186-113-128 128-186-113 128-113-186 113-186-128 113-128-186
'''





class Spectrum(object):
    """
    A class of the methods and variables that can analyze and describe a
    spectrum.
    Attributes:
        spectrum = the provided spectrum
        aa_mass = a dictionary of each amino acid's integer mass
    """
    def __init__(self, spectrum):
        self.spectrum = tuple(spectrum)
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


    def findMass(self, peptide):
        '''
        finds the mass of each peptide in the subpeptide list
        '''
        total = 0
        for p in peptide:
            total += self.aa_mass[p]
        return total


    def get_all_substrings(self,input_string):
        '''
        creates every possilbe substring of a string
        '''
        #stackoverflow
        length = len(input_string)
        return [input_string[i:j+1] for i in range(length) for j in range(i,length)]

    def checkConsistant(self, peptide):
        '''
        determines if the peptide is consistant with the provided spectrum.
        a peptide is consistant if its and all of its substrings' masses are
        within the spectrum.
        '''
        #if the peptide is greater than one, create all possible substrings
        #of the peptide and store it as a kmer
        if len(peptide) > 1:
            track = []
            kmers = []
            spec = list(self.spectrum)
            kmers = self.get_all_substrings(peptide)
            #iterate through all substrings of a peptide
            for k in kmers:
                mass = self.findMass(k)
                #if the mass of the substring is within the spectrum, track
                #it as True and remove it from the spectrum
                if mass in spec:
                    track.append(True)
                    spec.remove(mass)
                #if the mass of the substring is not within the spectrum, track
                #it as false
                else:
                    track.append(False)
            #if all of the substrings' masses are within the spectrum, return True
            if False not in track:
                return True
        else:
            track = []
            #repeats the same logic as above without making substrings
            if self.findMass(peptide) not in self.spectrum:
                track.append(False)
            else:
                track.append(True)
            if False not in track:
                return True
            else:
                return False


    def findCyclospectrum(self):
        '''
        Returns every amino acid string Peptide such that
        Cyclospectrum(Peptide) = Spectrum (if such a string exists)
        '''
        single_peps = self.makeAminoAcidMassDict()
        peps = self.makeAminoAcidMassDict()
        profiles = []
        #while the dictionary, peps, has values in it
        while bool(peps):
            #initiate a new dictionary
            new_peps = {}
            #iterate through the peps dictionary
            for p in peps:
                #iterate through the single peptide dictionary
                for s in single_peps:
                    #if the new peptide is consistant with the spectrum
                    if self.checkConsistant(p+s):
                        new_peps[p+s] = self.findMass(p+s)
            #iterate through the new peptides
            for n in new_peps:
                #if the new peptide's mass is equal to the maximum mass in the
                #spectrum
                if new_peps[n] == max(self.spectrum):
                    #append it to the profiles list
                    profiles.append(n)
            #set the dictionary peps equal to the new consistant peptides
            peps = new_peps
        #when the dictionary peps is empty, return the consistant profiles
        return profiles

    def printProfiles(self):
        '''
        prints the unique mass profiles of the amino acids in the
        Rosalind format
        '''
        #find the consistant profiles
        profiles = self.findCyclospectrum()
        profile_dict = {}
        profile_list = []
        formatted_list = []
        #iterate through each profile
        for p in profiles:
            #find the mass of each single peptide in the profile and store it
            #to a dictionary
            profile_dict[p] = [self.findMass(x) for x in p]
        #iterate through the dictionary of peptide profiles
        for p in profile_dict.values():
            #select only unique mass profiles and append it to a list
            if p not in profile_list:
                profile_list.append(p)
        #format the mass profiles to be separated by a "-"
        for p in profile_list:
            formatted_list.append("-".join('{0}'.format(n) for n in p))
        #print the formatted mass profiles as a space separated list
        print(" ".join(formatted_list))



def main():
    #the first list of data file is the spectrum
    s = sys.stdin.readline().strip().split()
    spectrum = []
    for i in s:
        spectrum.append(int(i))
    my_spectrum = Spectrum(spectrum)
    my_spectrum.printProfiles()


if __name__ == '__main__':
    main()
