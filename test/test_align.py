import pytest
from align import algs

@pytest.fixture

def some_relevant_data():
	return np.ones(10)

def test_fasta_io():
	test_object = algs.PairwiseAligner("scoring_matrices/BLOSUM62.mat", 10, 1)
	seq1 = test_object._read_sequence_file("sequences/prot-0004.fa")
	seq1_truth = "SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM"
	assert(seq1 == seq1_truth)

	seq2 = test_object._read_sequence_file("sequences/prot-4.fa")
	assert(seq2 == seq1_truth)
	
def test_scoring_matrix_io():
	test_object = algs.PairwiseAligner("scoring_matrices/BLOSUM62.mat", 10, 1)
	b62_dict = test_object.substitution_matrix
	b62_truth = {
		"A":{"A": 4, "R":-1, "N":-2, "D":-2, "C":0, "Q":-1, "E":-1, "G":0, "H":-2, "I":-1, "L":-1, "K":-1, "M":-1, "F":-2, "P":-1, "S":1, "T":0, "W":-3, "Y":-2, "V": 0, "B":-2, "Z":-1, "X":0, "*":-4},
		"R":{"A": -1, "R":5, "N":0, "D":-2, "C":-3, "Q":1, "E":0, "G":-2, "H":0, "I":-3, "L":-2, "K":2, "M":-1, "F":-3, "P":-2, "S":-1, "T":-1, "W":-3, "Y":-2, "V": -3, "B":-1, "Z":0, "X":-1, "*":-4},
		"N":{"A": -2, "R":0, "N":6, "D":1, "C":-3, "Q":0, "E":0, "G":0, "H":1, "I":-3, "L":-3, "K":0, "M":-2, "F":-3, "P":-2, "S":1, "T":0, "W":-4, "Y":-2, "V": -3, "B":3, "Z":0, "X":-1, "*":-4},
		"D":{"A": -2, "R":-2, "N":1, "D":6, "C":-3, "Q":0, "E":2, "G":-1, "H":-1, "I":-3, "L":-4, "K":-1, "M":-3, "F":-3, "P":-1, "S":0, "T":-1, "W":-4, "Y":-3, "V": -3, "B":4, "Z":1, "X":-1, "*":-4},
		"C":{"A": 0, "R":-3, "N":-3, "D":-3, "C":9, "Q":-3, "E":-4, "G":-3, "H":-3, "I":-1, "L":-1, "K":-3, "M":-1, "F":-2, "P":-3, "S":-1, "T":-1, "W":-2, "Y":-2, "V": -1, "B":-3, "Z":-3, "X":-2, "*":-4},
	}

	test_object = algs.PairwiseAligner("scoring_matrices/PAM250.mat", 10, 1)
	p250_dict = test_object.substitution_matrix
	p250_truth = {
		"A":{"A": 2, "R":-2, "N":0, "D":0, "C":-2, "Q":0, "E":0, "G":1, "H":-1, "I":-1, "L":-2, "K":-1, "M":-1, "F":-3, "P":1, "S":1, "T":1, "W":-6, "Y":-3, "V": 0, "B":0, "Z":0, "X":0, "*":-8},
		"R":{"A": -2, "R":6, "N":0, "D":-1, "C":-4, "Q":1, "E":-1, "G":-3, "H":2, "I":-2, "L":-3, "K":3, "M":0, "F":-4, "P":0, "S":0, "T":-1, "W":2, "Y":-4, "V": -2, "B":-1, "Z":0, "X":-1, "*":-8}
	}

	assert(b62_truth.items() <= b62_dict.items())
	assert(p250_truth.items() <= p250_dict.items())



def test_identical():
	assert True

def test_alignment_score():
	assert True


test_scoring_matrix_io()