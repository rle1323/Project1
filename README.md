# Project 1 - Sequence Alignment
## Due 01/27/2021

![BuildStatus](https://github.com/rle1323/Project1/workflows/HW1/badge.svg?event=push)

In this assignment, you will implement two classical alignment algorithms and then evaluate each algorithm’s performance with a range of parameters. There are two parts to this assignment and Part 2 requires completion of Part 1. We recommend reading through both Part 1 and Part 2 before beginning this assignment. 

* Part 1 - API and implementation
* Part 2 - Evaluating alignments

### main
Runs all code in align/\_\_main\_\_.py, useful for part 2
```
python -m align
```

### testing
Testing is as simple as running
```
python -m pytest test/*
```
from the root directory of this project.

### API Documentation

```
class PairwiseAligner
    This is the parent class for both pairwise alignment algorithms to be implemented in this assignment.
    It includes methods and attributes that are common to both algorithms.

    Parameters:
        sub_matrix_file::str
            Path to a substitution matrix file. This file is subsequently read into the substitution_matrix attribute
        gap_start_penalty::float
            This is the gap opening penalty for the affine gap implementation of local or global alignment.
        gap_extension_penalty::float
            This is the gap extension penalty for the affine gap implementation of local or global alignment
   
   
    method _read_substitution_matrix_file(self, substitution_matrix_filename):
        Reads in an amino acid substitution matrix from a .mat file, and stores substitution scores for 
        each amino acid pair in a 3d dictionary

        Arguments:
            substitution_matrix_filename::str
                Path to substitution matrix file
        
        Returns:
            substitution_dict::dict
                Dictionary with 24 primary keys corresponding to the 20 amino acids + 4 placeholders.
                Each primary key has 24 secondary keys. Indexing substitution_dict with a primary and secondary
                key will link to the substitution score for that given residue pair given the substitution matrix in use. 
  
  
    method _read_sequence_file(self, sequence_filename):
        Reads in a sequence file in FASTA format and returns it as a string.
        
        Arguments:
            sequence_filename::str
                Path to sequence file
        Returns:
            seq::str
                The protein sequence contained in the file
                
                
    method _initialize_score_matrix(self):
        Initializes an (m+1 x n+1) scoring/gap matrix full of zeros for use in either SW or NW, where m is the length of self.seqx and 
        n is the length of self.seqy
        
        Arguments:
            None
        
        Returns: 
            matrix::numpy array (int)
                Matrix of size (m x n) with only zeros as contents 
                
class SmithWaterman(PairwiseAligner):
    This class performs pairwise local alignment using the Smith-Waterman algorithm. It is a child class of PairwiseAligner, and therefore inherits
    the parameters necessary for initializing a SmithWaterman object.         
