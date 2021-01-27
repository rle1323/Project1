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

    method _forward_pass(self):
        This method performs the forward pass step of the Smith-Waterman algorithm, where all possible local alignments are scored and these 
        scores are held in a scoring matrix. A pointer matrix is also created to allow for the optimal local alignment to be found and constructed

        Arguments:
            None

        Returns:
            score_matrix::[float]
                A float array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the scores of all of the possible local alignments between the two sequences
            pointer_matrix_main::[int]
                An int array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the the traceback directions for each corresponding cell in score_matrix
                
                
    method _backwards_pass(self, pointer_matrix, i, j):
        This method performs traceback of any of the local alignments stored in score_matrix

        Arguments:
            pointer_matrix::[int]
                An int array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the the traceback directions for all possible alignment positions between the two sequences
            i::int
                starting position in seqx of the traceback being performed
            j::int
                starting position in seqy of the traceback being performed
        
        Returns:
            aligned_seqx::str
                The subsequence of seqx, with gaps, that is in the optimal local alignment between seqx and seqy
            aligned_seqy::str
                The subsequence of seqy, with gaps, that is in the optimal local alignment between seqx and seqy
                
                
    method align(self, seqx_file, seqy_file):
        Wrapper method that incorporates all of the steps for Smith-Waterman local alignment.

        Arguments:
            seqx_file::str
                Path to the fasta file containing the first sequence involved in the pairwise alignment
            seqy_file::str
                Path to the fasta file containing the second sequence involved in the pairwise alignment
        
        Returns:
            aligned_seqs::[str]
                A 1x2 list contaning both subsequences, with gaps if necessary, in the optimal local alignment
            alignment_score::float
                The alignment score of the optimal local alignment found between the two sequences      


class NeedlemanWunsch(PairwiseAligner):
    This class performs pairwise global alignment using the Needleman-Wunsch algorithm. It is a child class of PairwiseAligner, and therefore inherits
    the parameters necessary for initializing a NeedlemanWunsch object.

    method _forward_pass(self):
        This method performs the forward pass step of the Needleman-Wunsch algorithm, where all possible local alignments are scored and these 
        scores are held in a scoring matrix. A pointer matrix is also created to allow for the optimal global alignment to be constructed

        Arguments:
            None

        Returns:
            score_matrix::[float]
                A float array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the scores of all of the possible global alignments between the two sequences
            pointer_matrix_main::[int]
                An int array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the the traceback directions for each corresponding cell in score_matrix
                
                
    method _backwards_pass(self, pointer_matrix, i, j):
        This method performs traceback of the optimal global alignment stored in score_matrix

        Arguments:
            pointer_matrix::[int]
                An int array of size (m+1)x(n+1), where m and n are the lengths of the two sequences being aligned respectively. This matrix
                holds the the traceback directions for all possible alignment positions between the two sequences
            i::int
                starting position in seqx of the traceback being performed. For Needleman-Wunsch, this is always equal to m
            j::int
                starting position in seqy of the traceback being performed. For Needleman-Wunsch, this is always equal to n
        Returns:
            aligned_seqx::str
                The sequence of seqx, with gaps, that is in the optimal global alignment between seqx and seqy
            aligned_seqy::str
                The sequence of seqy, with gaps, that is in the optimal global alignment between seqx and seqy
          
          
    method align(self, seqx_file, seqy_file):
        Wrapper method that incorporates all of the steps for Needleman-Wunsch global alignment.

        Arguments:
            seqx_file::str
                Path to the fasta file containing the first sequence involved in the pairwise alignment
            seqy_file::str
                Path to the fasta file containing the second sequence involved in the pairwise alignment
        
        Returns:
            aligned_seqs::[str]
                A 1x2 list contaning both sequences, with gaps if necessary, in the optimal global alignment
            alignment_score::float
<<<<<<< HEAD
                The  score of the optimal global alignment found between the two sequences
=======
                The  score of the optimal global alignment found between the two sequences
>>>>>>> a56d673602b7436975e4f8ec33d25169ad565411
