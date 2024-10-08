Steps taken in build a CNV phylogenetic tree algorithm.

We base ourselves from the phylogenetic tree algorithm from sequence data with bayesian inference. We first develop an algorithm for sequence data with Metropolis-Hastings MCMC, because of simplicity. We then modify the data representation, likelihood function, distance metric, tree proposal metric. 

Interpretation of CNV Values

0.0: Complete deletion of the gene in both DNA strands
0.5: Deletion in one DNA strand
1.0: Baseline (normal copy number)
1.5: Amplification in one DNA strand
2.0: Amplification in both DNA strands
2.5: Three copies (amplification in both strands + additional copy in one strand)
3.0: Four copies (amplification in both strands + additional copy in both strands)

Example Interpretations

Species A [1.0, 1.5, 1.0, 0.5, 1.0, 1.0, 2.0]:

Normal, Amplified in one strand, Normal, Deleted in one strand, Normal, Normal, Amplified in both strands


Species D [1.5, 1.5, 0.5, 1.0, 1.5, 2.0, 1.0]:

Amplified in one strand, Amplified in one strand, Deleted in one strand, Normal, Amplified in one strand, Amplified in both strands, Normal


Hypothetical example with 0.0 [0.0, 1.0, 2.0, 1.5]:

Complete deletion, Normal, Amplified in both strands, Amplified in one strand

For the single cell example.

But now we want to change this to spot data, so they should be the average of single cell in a spot.