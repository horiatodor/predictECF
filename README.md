# predictECF

Author: Horia Todor [horia.todor@gmail.com]

## How to use this code

There are two main functions, `predict35w` and  `predict10w`, both have the same input format:

`predict35w(sequence, dnamatrix, proteinmatrix, row_weights=NA)`

`sequence` is a vector containing the protein sequence of the ECF sigma factor whose promoter specificity is being predicted. This should be a snippet of an alignment that matches supplementary table 1. 

`dnamatrix` and `proteinmatrix` are matched matrices of the protein and dna alignments. The package includes these as lazy data. For 'predict35w' these should be respectively `motif35` and `domain4`.  For 'predict10w' these should be respectively `motif10` and `domain2`.  

`row_weights` is a vector of the same length as `dnamatrix` and contains the relative weights of each pair of protein/dna sequences. The package includes 'applied_weights_3' as lazy data: this should be used with `motif35` and `domain4` or `motif10` and `domain2`.  

Output is a 4x7 matrix containing the log likelyhood of each DNA base at each position based on the alignments input.   
