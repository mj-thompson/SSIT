# SSIT

Stringent, Sequencing and Inversion Technique

SSIT was a final project that I worked on as an undergraduate in a bioinformatics algorithms course. From a very basic alignment algorithm, I changed a bit of the functionality to include alignment by hashing, SNP-detection, and inversion detection.

The algoirthm essentially uses a low threshold of mismatches to call SNPs, and record areas of heavy mutation (STR, CNV, indels, inversions). For the sake of the final, I only searched for inversions. While the algorithm has a high number of false positives (the filtering script attempts to mitigate this issue) and can really be improved, it provides great insight into regions of mutation when comparing an individual to a reference genome.
