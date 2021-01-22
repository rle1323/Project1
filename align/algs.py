
class PairwiseAligner:
	"""
	This is the parent class for both pairwise alignment algorithms to be implemented in this assignment.
	It includes methods and attributes that are common to both algorithms.
	"""
	def __init__(self, seqx_file, seqy_file, sub_matrix_file, gap_start_penalty, gap_extension_penalty):
		self.seqx
		pass
	
	def _read_substitution_matrix_file(substitution_matrix_filename):
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
		with open(substitution_matrix_file) as f:
			residue_counter = 0
			for line in f.readlines():
				line_items = line.strip().split()
				# handle comment lines by simply ignoring them (return to loop)
				if line_items[0] == "#" or line_items[0] == ".":
					continue
				# if we have the list of residues, save them
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

	def _read_sequence_file(sequence_filename):
		"""
		Reads in a sequence file in FASTA format

		Arguments:
			sequence_filename::str
				Path to sequence file
		Returns:
			seq::str
				The protein sequence contained in the file
		"""
		with open(sequence_filename) as f:
			seq = ""
			for line in f.readlines():
				line_item = line.strip()
				if line_item.startswith(">") is False: 
					seq += line_item
		return seq


class SmithWaterman(PairwiseAligner):
	pass

class NeedlemanWunsch(PairwiseAligner):
	pass