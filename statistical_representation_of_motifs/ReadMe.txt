To determine which motif in a sequence is the most statistically underrepresented and therefore determine significant information about the sequence such as the exact role of a non-coding RNA in gene expression(Cserzo, 2010), I analyzed kmers of varying lengths from select parts of a specie's genome. 

To determine the most statistically underrepresented kmer, I calculated the actual number of occurrences, the expected number of occurrences, and the Z-score of each kmer. From the Z-score values, one can easily identify which kmer statistically occurs less often than expected. The most statistically underrepresented motif is easily identified as the program sorts the output by the length of the kmer, then by the Z-score. 


References:
Cserzo, Miklos et al. “Relating underrepresented genomic DNA patterns and tiRNAs: the rule behind the observation and beyond” Biology direct vol. 5 56. 22 Sep. 2010, doi:10.1186/1745-6150-5-56