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
