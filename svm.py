
from xlrd import open_workbook
from sklearn import svm
from math import pow
from collections import defaultdict 
import hashlib
import itertools

#http://en.wikipedia.org/wiki/Amino_acid#Classification
"""
Aliphatic:	Glycine, Alanine, Valine, Leucine, Isoleucine
Hydroxyl or Sulfur/Selenium-containing:	Serine, Cysteine, Selenocysteine, Threonine, Methionine
Cyclic:	Proline
Aromatic:	Phenylalanine, Tyrosine, Tryptophan
Basic:	Histidine, Lysine, Arginine
Acidic and their Amide:	Aspartate, Glutamate, Asparagine, Glutamine
"""
#Map using above grouping of AA's
AA1 = {"G":1, "A":1, "V":1, "L":1, "I":1,
		"S":2, "C":2, "U":2, "T":2, "M":2,
		"P":3,
		"F":4, "Y":4, "W":4,
		"H":5, "K":5, "R":5,
		"D":6, "E":6, "N":6, "Q":6}

#http://en.wikipedia.org/wiki/Proteinogenic_amino_acid#Non-specific_abbreviations
AA2 = {'R':1, 'H':1, 'K':1,
		'D':2, 'E':2,
		'S':3, 'T':3, 'N':3, 'Q':3,
		'C':4, 'U':4, 'G':4, 'P':4,
		'A':5, 'I':5, 'L':5, 'M':5, 'F':5, 'W':5, 'Y':5, 'V':5}

#http://www.ann.com.au/MedSci/amino.htm (for AA3, AA4, AA5)
# nonpolar, polar, acidic(polar), basic(polar)
AA3 = {'G':1, 'A':1, 'V':1, 'L':1, 'I':1, 'P':1, 'M':1, 'F':1, 'W':1,
	'S':2, 'T':2, 'N':2, 'Q':2, 'C':2, 'Y':2,
	'D':3, 'E':3,
	'K':4, 'R':4, 'H':4}
#structure
AA4 = {'G':1, 'A':1, 
		'V':2, 'L':2, 'I':2, 
		'P':3, 'F':3, 
		'Y':4, 'W':4, 
		'M':5, 
		'S':6, 'T':6, 
		'C':7,
		'Q':8, 'N':8,
		'D':9, 'E':9,
		'K':10, 'R':10, 'H':10
		}
AA5 = {'G':1, 'A':1, 'V':1, 'L':1, 'I':1,
		'C':2, 'M':2,
		'F':3, 'Y':3, 'W':3,
		'S':4, 'T':4, 'N':4, 'Q':4,
		'E':5, 'D':5,
		'K':6, 'R':6, 'H':6,
		'P':7 
		}	

"""
G - Glycine (Gly)
P - Proline (Pro)
A - Alanine (Ala)
V - Valine (Val)
L - Leucine (Leu)
I - Isoleucine (Ile)
M - Methionine (Met)
C - Cysteine (Cys)
F - Phenylalanine (Phe)
Y - Tyrosine (Tyr)
W - Tryptophan (Trp)
H - Histidine (His)
K - Lysine (Lys)
R - Arginine (Arg)
Q - Glutamine (Gln)
N - Asparagine (Asn)
E - Glutamic Acid (Glu)
D - Aspartic Acid (Asp)
S - Serine (Ser)
T - Threonine (Thr)
"""
#One letter codes for AA's
AA = ['G', 'P', 'A', 'V', 'L', 'I', 'M', 'C',
		'F', 'Y', 'W', 'H', 'K', 'R', 'Q', 'N',
		'E', 'D', 'S', 'T']

#Returns list of all possible n_grams
def all_n_grams(n):
	#each elt. in list is a tuple of AA's
	grams = list(itertools.product(AA, repeat=n))
	
	#convert tuples to strings
	for index, AA_tuple in enumerate(grams):
		gram = ""
		for a in AA_tuple:
			gram += a
		grams[index] = gram
	return grams

#Returns vector counting how many of each possible n_gram occured
#Length of vector is length of all_n_grams
#occurences[i] = # times all_n_grams[i] occured in seq
def n_gram_counts(seq, n, all_n_grams):
	#n_grams = all_n_grams(n)
	occurences = [0]*len(all_n_grams)
	for i in range(0, len(seq) - n):
		gram = ''
		for j in range(i, i+n):
			gram += seq[j]
		index = all_n_grams.index(gram)
		occurences[index] += 1
	return occurences

#Count occurences of each Amino Acid
def AA_counts(seq):
	AA = {'G':0, 'P': 1, 'A': 2, 'V':3, 'L':4,
		'I':5, 'M':6, 'C':7, 'F':8, 'Y':9,
		'W':10, 'H':11, 'K':12, 'R':13, 'Q':14,
		'N':15,'E':16, 'D':17, 'S':18, 'T':19}
	counts = [0] * 20
	for a in seq:
		counts[AA[a]] += 1
	return counts

#Count occurences of each class of AA
def AAn_counts(seq, AAn):
	counts = [0] * len(set(AAn.values()))
	for a in seq:
		#Grouping values consecutive indexed by 1
		counts[AAn[a] - 1] += 1
	return counts 


max_seq_len = 465 
#Map sequence to groupings
#Len(feature vector) = len(seq)
#v[i] = which group seq[i] is in according to AAn
def map_seq(seq, AAn):
	mapped_features = []
	for a in seq: #map AA's using 1st AA grouping
		mapped_features.append(AAn[a])
	mapped_features += [0 for _ in range(max_seq_len-len(seq))] # pad with zero's (seqs have diff lengths)
	return mapped_features

#Map for different substrates (classes)
subs_class_dict = {'dhpg':0,'horn':1, 'pip':2, 'bht':3, 'dab':4,'dhb':5,
		'Orn':6, 'dht':7, 'hpg':8, 'A':9, 'C':10, 'E':11, 'D':12, 'G':13,
		'F':14, 'I':15, 'K':16, 'L':17, 'N':18, 'Q':19, 'P':20, 'S':21,
		'R':22,'T':23, 'W':24, 'V':25, 'Y':26, 'orn':27,
		'beta-ala':28, 'ORN':29, 'hyv-d':30, 'aad':31}

#Open Workbook, return tuple (vector of sequences, vector of substrate classes)
def getData():
	seqs=[]
	subs=[]

	book = open_workbook('Adomain_Substrate.xls')
	worksheet = book.sheet_by_name('Adomain_Substrate')
	num_rows = worksheet.nrows
	substrate_cell = 1
	seq_cell = 2
	for i in range (1, num_rows):
		seq = worksheet.cell_value(i, seq_cell)
		substrate = worksheet.cell_value(i, substrate_cell)
		seqs.append(seq)
		subs.append(subs_class_dict[substrate])
	return (seqs, subs)

#Return feature vector for a sequence
#Takes variable keyword arguments
#feature ="ngram", "AAcounts", "AAncounts", or "mapseq"
#extra parameter n depending on given feature
def getFeaturesFromSeq(seq, **kwargs):
	features = []
	f = kwargs["feature"]
	AAn_dict = {1: AA1, 2: AA2, 3:AA3, 4:AA4, 5:AA5}
	if f == "ngram":
		n = kwargs["n"]
		grams = all_n_grams(n)
		features = n_gram_counts(seq, n, grams)
	elif f == "AAcounts":
		features = AA_counts(seq)
	elif f == "AAncounts":
		n = kwargs["n"]
		features = AAn_counts(seq, AAn_dict[n])
	elif f == "mapseq":
		n = kwargs["n"]
		features = map_seq(seq, AAn_dict[n])
	return features

#Adds features to pre-existing X using getFeaturesFromSeq(seq, **kwargs)
def addFeatures(seqs, X, **kwargs):
	for i, seq in enumerate(seqs):
		X[i] += getFeaturesFromSeq(seq, **kwargs)
	return X


"""  EXAMPLE

data = getData()
seqs = data[0]
X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
Y = data[1]

X = addFeatures(seqs, X, feature="ngram", n=2)
X = addFeatures(seqs, X, feature="mapseq", n=1)

X_train = X[0:1100]
Y_train = Y[0:1100]
X_test = X[1100:]
Y_test = Y[1100:]

clf= svm.SVC(kernel='linear').fit(X_train, Y_train)
print clf.score(X_test, Y_test)

"""