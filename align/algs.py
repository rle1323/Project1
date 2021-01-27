import numpy as np
import math

class PairwiseAligner:
    """
    This is the parent class for both pairwise alignment algorithms to be implemented in this assignment.
    It includes methods and attributes that are common to both algorithms.

    Parameters:
        sub_matrix_file::str
            Path to a substitution matrix file. This file is subsequently read into the substitution_matrix attribute
        gap_start_penalty::float
            This is the gap opening penalty for the affine gap implementation of local or global alignment.
        gap_extension_penalty::float
            This is the gap extension penalty for the affine gap implementation of local or global alignment
    """

    def __init__(self, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
        self.seqx = None
        self.seqy = None
        self.substitution_matrix = self._read_substitution_matrix_file(sub_matrix_file)
        self.gap_start_penalty = gap_start_penalty
        self.gap_extension_penalty = gap_extension_penalty
    
    def _read_substitution_matrix_file(self, substitution_matrix_filename):
        """
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
        """
        with open(substitution_matrix_filename) as f:
            residue_counter = 0
            for line in f.readlines():
                # remove unnecessary space (including "\n") and split the contents of the line into a list
                line_items = line.strip().split()
                # handle comment lines by simply ignoring them (return to loop)
                if line_items[0] == "#" or line_items[0] == ".":
                    continue
                # if we have the list of residues, save them and initialize a dictionary with each residue as a primary key
                elif line_items[0] == "A":
                    residues = line_items
                    substitution_dict = {residue : None for residue in residues}
                # otherwise, line_items is a list of scores for a particular residue with the other residues
                else:
                    # get the residue for which we are inputting the rest of the scores
                    resid = residues[residue_counter]
                    residue_counter += 1
                    # fill in that dictionary slot with the corresponding scores
                    substitution_dict[resid] = {residues[i] : int(line_items[i]) for i, residue in enumerate(line_items)}
        return substitution_dict

    def _read_sequence_file(self, sequence_filename):
        """
        Reads in a sequence file in FASTA format

        Arguments:
            sequence_filename::str
                Path to sequence file
        Returns:
            seq::str
                The protein sequence contained in the file
        """
        # open file, read in line by line. If line does not start with ">" (fasta header marker), add it to the sequence
        with open(sequence_filename) as f:
            seq = ""
            for line in f.readlines():
                line_item = line.strip()
                if line_item.startswith(">") is False: 
                    seq += line_item
        seq = seq.upper()
        return seq

    def _initialize_score_matrix(self):
        """
        Initializes an (m+1 x n+1) scoring/gap matrix full of zeros for use in either SW or NW, where m is the length of self.seqx and 
        n is the length of self.seqy
        
        Arguments:
            None
        
        Returns: 
            matrix::numpy array (int)
                Matrix of size (m x n) with only zeros as contents 
        """
        matrix = np.zeros( (len(self.seqx) + 1, len(self.seqy) + 1) )
        return matrix 

    def _forward_pass(self):
        """
        Generic parent version of the method that fills in the scoring matrix for either local or global alignment. Class-specific implementations included in
        class definitions.
        """
        pass

    def _backwards_pass(self):
        """
        Generic parent version of the method that follows pointer matrix to reconstruct alignment for either local or global alignment. Class-specific implementations included in
        class definitions.
        """
        pass
    
    def align(self):
        """
        Generic parent version of the wrapper method that performs entire alignment process. Class-specific implementations of this method 
        are largely the same, with minor changes in how alignment scores are found, and where to start alignment traceback
        """
        pass

class SmithWaterman(PairwiseAligner):
    """
    This class performs pairwise local alignment using the Smith-Waterman algorithm. It is a child class of PairwiseAligner, and therefore inherits
    the parameters necessary for initializing a SmithWaterman object.
    """

    def __init__(self, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
        super().__init__(sub_matrix_file, gap_start_penalty, gap_extension_penalty)

    def _forward_pass(self):
        """
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
        """
        # initialize matrices for storing scores 
        score_matrix = self._initialize_score_matrix()
        gap_matrix_x = self._initialize_score_matrix()
        gap_matrix_y = self._initialize_score_matrix() 
        

        # for SW affine gap alignment scheme, need to make first row and first column of gap matrices = -infinite
        gap_matrix_x[:,0] = float("-inf")
        gap_matrix_y[:,0] = float("-inf")
        gap_matrix_x[0,:] = float("-inf")
        gap_matrix_y[0,:] = float("-inf")

        # initialize pointer matrix
        pointer_matrix_main =  self._initialize_score_matrix()
        # fill out score matrices and traceback matrix (forward)
        for i in range(score_matrix.shape[0]):
            if i == 0: continue # leave first column with zeros
            for j in range(score_matrix.shape[1]):
                if j == 0: continue # leave first row with zeros
                # since dims of the scoring matrix are 1 greater than the sequence lengths, correct for that and load the respective residues
                resid_x = self.seqx[i-1]
                resid_y = self.seqy[j-1]
                # retrieve the substitution/match score of the current residues in the sequences
                substitution_score = self.substitution_matrix[resid_x][resid_y]
                # if regular forward step, update the reccurence rule matrices

                # update gap matrices
                gap_matrix_x[i,j] = max(gap_matrix_x[i,j-1] - self.gap_extension_penalty, # extend gap in x
                                        score_matrix[i,j-1] - self.gap_start_penalty) #open gap in x
                
                gap_matrix_y[i,j] = max(gap_matrix_y[i-1,j] - self.gap_extension_penalty, # extend gap in y
                                        score_matrix[i-1,j] - self.gap_start_penalty) # open gap in y

                # update main scoring matrix
                score_matrix[i,j] = max(score_matrix[i-1,j-1] + substitution_score,
                                        gap_matrix_x[i,j],
                                        gap_matrix_y[i,j],
                                        0)			

                # update main pointer matrix, where 3 indicates to stop the traceback
                if score_matrix[i,j] == (score_matrix[i-1,j-1] + substitution_score): pointer_matrix_main[i,j] = 0
                elif score_matrix[i,j] == gap_matrix_x[i,j]: pointer_matrix_main[i,j] = 1
                elif score_matrix[i,j] == gap_matrix_y[i,j]: pointer_matrix_main[i,j] = 2
                else: pointer_matrix_main[i,j] = 3 # at this position, main score matrix == 0, this marks to stop the alignment
        
        return score_matrix, pointer_matrix_main
    
    def _backwards_pass(self, pointer_matrix, i, j):
        """
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
        """
        aligned_seqx = ""
        aligned_seqy = ""
        while i > 0 and j > 0 and pointer_matrix[i,j] != 3:
            current_val = pointer_matrix[i,j]
            if current_val == 0 : # this means that the previous move in the alignment came from dai
                aligned_seqx = self.seqx[i-1] + aligned_seqx
                aligned_seqy = self.seqy[j-1] + aligned_seqy
                i -= 1
                j -= 1
            elif current_val == 1: # this means we are in the x axis scoring matrix
                aligned_seqx = "-" + aligned_seqx
                aligned_seqy = self.seqy[j-1] + aligned_seqy
                j -= 1 # since we have a gap in x, only go to the next residue for y
            elif current_val == 2: # this means that we are in the y-axis pointer matrix
                aligned_seqy = "-" + aligned_seqy
                aligned_seqx = self.seqx[i-1] + aligned_seqx
                i -= 1 # since we have a gap in y, only go to the next residue for x
        return aligned_seqx, aligned_seqy
    
    def align(self, seqx_file, seqy_file):
        """
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
        """
        self.seqx = self._read_sequence_file(seqx_file)
        self.seqy = self._read_sequence_file(seqy_file)
        # fill out score matrices and traceback matrix (forward)
        score_matrix, pointer_matrix_main = self._forward_pass()

        max_i, max_j = np.unravel_index(score_matrix.argmax(), score_matrix.shape)
        alignment_score = np.amax(score_matrix)
        # traceback (backwards pass)	
        # get index of maximum score in score_matrix, to start the traceback at that index
        aligned_seqx, aligned_seqy = self._backwards_pass(pointer_matrix_main, max_i, max_j)
        aligned_seqs = [aligned_seqx, aligned_seqy]
        return aligned_seqs, alignment_score
       

class NeedlemanWunsch(PairwiseAligner):
    """
    This class performs pairwise global alignment using the Needleman-Wunsch algorithm. It is a child class of PairwiseAligner, and therefore inherits
    the parameters necessary for initializing a NeedlemanWunsch object.
    """

    def __init__(self, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
        super().__init__(sub_matrix_file, gap_start_penalty, gap_extension_penalty)

    def _forward_pass(self):
        """
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
        """
        # initialize matrices for storing scores 
        score_matrix = self._initialize_score_matrix()
        gap_matrix_x = self._initialize_score_matrix()
        gap_matrix_y = self._initialize_score_matrix() 

        # for NW affine gap alignment scheme, need to make first row and first column of gap matrices = -i*GE - GO
        for i in range(score_matrix.shape[0]):
            score_matrix[i,0] = i*(-self.gap_extension_penalty) - self.gap_start_penalty
            gap_matrix_x[i,0] = i*(-self.gap_extension_penalty) - self.gap_start_penalty
        for j in range(score_matrix.shape[1]):
            score_matrix[0,j] = j*(-self.gap_extension_penalty) - self.gap_start_penalty
            gap_matrix_y[0,j] = j*(-self.gap_extension_penalty) - self.gap_start_penalty
        
        score_matrix[0,0] = 0

        # initialize pointer matrix
        pointer_matrix_main =  self._initialize_score_matrix()
        # fill out score matrices and traceback matrix (forward)
        for i in range(score_matrix.shape[0]):
            if i == 0: continue # leave first column initial values
            for j in range(score_matrix.shape[1]):
                if j == 0: continue # leave first row with initial values
                # since dims of the scoring matrix are 1 greater than the sequence lengths, correct for that and load the respective residues
                resid_x = self.seqx[i-1]
                resid_y = self.seqy[j-1]
                # retrieve the substitution/match score of the current residues in the sequences
                substitution_score = self.substitution_matrix[resid_x][resid_y]
                # if regular forward step, update the reccurence rule matrices
                # update gap matrices
                gap_matrix_x[i,j] = max(gap_matrix_x[i,j-1] - self.gap_extension_penalty, # extend gap in x
                                        score_matrix[i,j-1] - self.gap_start_penalty) #open gap in x
                
                gap_matrix_y[i,j] = max(gap_matrix_y[i-1,j] - self.gap_extension_penalty, # extend gap in y
                                        score_matrix[i-1,j] - self.gap_start_penalty) # open gap in y

                # update main scoring matrix
                score_matrix[i,j] = max(score_matrix[i-1,j-1] + substitution_score,
                                        gap_matrix_x[i,j],
                                        gap_matrix_y[i,j])
                
                # update pointer matrix
                pointer_matrix_main[i,j] = np.argmax((score_matrix[i-1,j-1] + substitution_score,
                                                      gap_matrix_x[i,j],
                                                      gap_matrix_y[i,j]))
        return score_matrix, pointer_matrix_main
    
    def _backwards_pass(self, pointer_matrix, i, j):
        """
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
        """
        aligned_seqx = ""
        aligned_seqy = ""
        while i > 0 or j > 0:
            current_val = pointer_matrix[i,j]
            # cut off the traceback. Fill in the rest of the unfinished sequence with its preceding portion an put in gaps for the other seq
            if i <= 0:
                aligned_seqy = self.seqy[:j] + aligned_seqy
                aligned_seqx = "-"*len(self.seqy[:j]) + aligned_seqx
                break
            # cut off the traceback. Fill in the rest of the unfinished sequence with its preceding portion an put in gaps for the other seq
            elif j <= 0:
                aligned_seqx = self.seqx[:i] + aligned_seqx
                aligned_seqy = "-"*len(self.seqx[:i]) + aligned_seqy
                break
            # otherwise, continue as normal
            if current_val == 0 : # this means that we are in the main (diag) scoring matrix
                aligned_seqx = self.seqx[i-1] + aligned_seqx
                aligned_seqy = self.seqy[j-1] + aligned_seqy
                i -= 1
                j -= 1
            elif current_val == 1: # this means we are in the x axis scoring matrix
                aligned_seqx = "-" + aligned_seqx
                aligned_seqy = self.seqy[j-1] + aligned_seqy
                j -= 1 # since we have a gap in x, only go to the next residue for y
            elif current_val == 2: 
                aligned_seqy = "-" + aligned_seqy
                aligned_seqx = self.seqx[i-1] + aligned_seqx
                i -= 1 # since we have a gap in y, only go to the next residue for x
        return aligned_seqx, aligned_seqy

    
    def align(self, seqx_file, seqy_file):
        """
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
                The  score of the optimal global alignment found between the two sequences
        """
        # read in the two sequence files
        self.seqx = self._read_sequence_file(seqx_file)
        self.seqy = self._read_sequence_file(seqy_file)
        # fill out score matrices and traceback matrix (forward)
        score_matrix, pointer_matrix_main = self._forward_pass()
        alignment_score = score_matrix[-1,-1]
        # traceback (backwards pass)	
        # get index of maximum score in score_matrix, to start the traceback at that index
        aligned_seqx, aligned_seqy = self._backwards_pass(pointer_matrix_main, len(self.seqx), len(self.seqy))
        aligned_seqs = [aligned_seqx, aligned_seqy]
        return aligned_seqs, alignment_score