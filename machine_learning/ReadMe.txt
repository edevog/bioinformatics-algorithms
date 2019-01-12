Hidden Markov models (HMMs) have been used to identify protein coding genes in organisms' DNA, such as E. coli (Krogh, 1994). This program is an exercise to understand this method by first applying it to a simpler version of the problem.

HMM_path and HMM_probability use HMMs to find the most probable hidden path for a given sequence and the probability of the HMM emitting a specific string respectively. 

To implement more accurate HMMs where the transition and emission probabilities are unknown, I implemented two different methods of learning---Viterbi and Baum-Welch---in order to determine the transition and emission probabilities.



References:
Krogh, A et al. “A hidden Markov model that finds genes in E. coli DNA” Nucleic acids research vol. 22,22 (1994): 4768-78.