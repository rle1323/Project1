import numpy as np
import math

class PairwiseAligner:
	"""
	This is the parent class for both pairwise alignment algorithms to be implemented in this assignment.
	It includes methods and attributes that are common to both algorithms.
	"""
	def __init__(self, seqx_file, seqy_file, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
		self.seqx = None
		self.seqy = None
		self.substitution_matrix = None
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
		return seq

	def _initialize_score_matrix(self):
		"""
		Initializes an (m x n) scoring/gap matrix full of zeros for use in either SW or NW, where m is the length of self.seqx and 
		n is the length of self.seqy
		
		Arguments:
			None
		
		Returns: 
			matrix::numpy array (int)
				Matrix of size (m x n) with only zeros as contents 
		"""
		matrix = np.zeros( (len(self.seqx), len(self.seqy)) )
		return matrix 

class SmithWaterman(PairwiseAligner):
	def __init__(self, seqx_file, seqy_file, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
		super().__init__(seqx_file, seqy_file, sub_matrix_file, gap_start_penalty, gap_extension_penalty)
		self.seqx = self._read_sequence_file(seqx_file)
		self.seqy = self._read_sequence_file(seqy_file)
		self.substitution_matrix = self._read_substitution_matrix_file(sub_matrix_file)
	
	def align(self):

		# initialize matrices for storing scores 
		score_matrix = self._initialize_score_matrix()
		gap_matrix_x = self._initialize_score_matrix()
		gap_matrix_y = self._initialize_score_matrix() 

		# initialize pointer matrices
		score_pointer_matrix =  self._initialize_score_matrix()
		gap_pointer_matrix_x = self._initialize_score_matrix()
		gap_pointer_matrix_y = self._initialize_score_matrix() 

		# fill out score matrices and traceback matrix
		for i, resid_x in enumerate(self.seqx):
			if i == 0: continue # leave first column with zeros
			for j, resid_y in enumerate(self.seqy):
				if j == 0: continue # leave first row with zeros
				# update gap matrices
				gap_matrix_x[i,j] = max(gap_matrix_x[i-1,j] - self.gap_extension_penalty, # continue gap in seqy
										score_matrix[i-1,j] - (self.gap_start_penalty + self.gap_extension_penalty),
										0) #start gap in w
				
				gap_matrix_y[i,j] = max(gap_matrix_y[i,j-1] - self.gap_extension_penalty,
										score_matrix[i,j-1] - (self.gap_start_penalty + self.gap_extension_penalty),
										0)
				# update main scoring matrix
				score_matrix[i,j] = max(score_matrix[i-1, j-1] + self.substitution_matrix[resid_x][resid_y],
										gap_matrix_x[i,j],
										gap_matrix_y[i,j],
										0)
															 
		print(np.amax(score_matrix), np.amax(gap_matrix_x), np.amax(gap_matrix_y))


				
	
class NeedlemanWunsch(PairwiseAligner):
	pass

sw_test = SmithWaterman(seqx_file = "../sequences/prot-0004.fa",
						seqy_file = "../sequences/prot-0014.fa",
						sub_matrix_file = "../scoring_matrices/BLOSUM62.mat",
						gap_start_penalty = 10,
						gap_extension_penalty = 1)

print(sw_test.seqx, "\n")
print(sw_test.seqy)
sw_test.align()
