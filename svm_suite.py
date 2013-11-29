from svm import *
from xlrd import open_workbook
from xlwt import *
from scrape_uniprot import *


#example
#opens workbooks
book = open_workbook('Adomain_Substrate.xls')
worksheet = book.sheet_by_name('Adomain_Substrate')

book_write = Workbook()
data = getData()

def poly_dimension_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]

	X = addFeatures(seqs, X, feature="ngram", n=2)
	X = addFeatures(seqs, X, feature="mapseq", n=1)
	worksheet_write = book_write.add_sheet("Linear Degrees")
	write(0, 0, "Linear Degrees", worksheet_write)
	for i in range(1, 10):
		clf = clf= svm.SVC(kernel='poly', degree=i)
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', k=3)
		write(i, 0, val, worksheet_write)

def experiment_with_k():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]

	X = addFeatures(seqs, X, feature="ngram", n=2)
	worksheet_write = book_write.add_sheet("Kernels")
	write(0, 1, "K-Varies", worksheet_write)

	for i in range(2, 10):
		#clf = clf= svm.SVC(kernel='poly', degree=)
		clf = clf= svm.SVC(kernel='linear')
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', k=i)
		write(i-1, 1, val, worksheet_write)


# Feature experiment
# Linear kernel, k=4

def feature_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]
	worksheet_write = book_write.add_sheet("features")
	#write(0,0,"Features", worksheet_write)
	clf = clf= svm.SVC(kernel='linear')
	row = 0

	"""
	#AAcounts
	X = addFeatures(seqs, X, feature="AAcounts")
	val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
	write(row, 1, val, worksheet_write)
	write(row, 0, "AAcounts", worksheet_write)
	row+=1

	#AAncounts
	for i in range(1,6):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="AAncounts", n=i) #n=1-5
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(row, 0, "AAncounts, n=" + str(i), worksheet_write)
		write(row, 1, val, worksheet_write)
		row+=1
	#mapseq
	for i in range(1,6):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="mapseq", n=i) #n=1-5
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(row, 0, "mapseq, n=" + str(i), worksheet_write)
		write(row, 1, val, worksheet_write)
		row+=1

	"""
	#ngrams
	for i in range(2,3):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="ngram", n=i) 
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(row, 0, "ngrams, n=" + str(i), worksheet_write)
		write(row, 1, val, worksheet_write)
		row+=1
	"""

	#AAn_ngrams
	for i in range(1, 6):
		for j in range (3, 5):
			X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
			X = addFeatures(seqs, X, feature="AAn_ngram", AAn=i, n=j)
			val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
			write(row, 0, "AAn_ngram, AAn=" + str(i) + ", n=" + str(j)
				, worksheet_write)
			write(row, 1, val, worksheet_write)
			row+=1
	
	#AA_distances
	for i in range(21, 31):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="AA_distances", n=i)
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(row, 0, "AA_distances, n=" + str(i), worksheet_write)
		write(row, 1, val, worksheet_write)
		row+=1

	#AAn_distances
	for i in range(1, 11):
		for j in range(1,6):
			X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
			X = addFeatures(seqs, X, feature="AAn_distances", n=i, AAn=j)
			val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
			write(row, 0, "AAn_distances, n=" + str(i) +", AAn=" + str(j),
				worksheet_write)
			write(row, 1, val, worksheet_write)
			row+=1
	"""

def rbf_gridsearch_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]

	worksheet_write = book_write.add_sheet("rbf_gridsearch")

	c_values = [pow(2,i) for i in range(-5,16)[::2]]
	gamma_values = [pow(2,i) for i in range(-15,4)[::2]]

	#write first row c headers to excel file
	col = 1
	for col, c in enumerate(c_values):
		write(0, col, "C=" + str(c), worksheet_write)
	
	#gridsearch
	row = 1
	for gamma in gamma_values:
		col = 0
		#write left column header gamma
		write(row, col, "g=" + str(gamma), worksheet_write)
		for c in c_values:
			col += 1
			clf = clf= svm.SVC(kernel='rbf', C=c, gamma=gamma)
			X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
			X = addFeatures(seqs, X, feature="ngram", n=2)
			val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
			write(row, col, val, worksheet_write)
		row+=1

from sklearn.metrics import *
from sklearn import svm
import numpy as np

def confusion_matrix_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	X = addFeatures(seqs, X, feature="ngram", n=2)
	Y = data[1]
	clf = clf= svm.SVC(kernel='linear')

	folds = 4
	skf = StratifiedKFold(np.array(Y), folds)
	i = 0
	for train, test in skf:
		clf.fit(np.array(X)[train], np.array(Y)[train])
		Y_true = np.array(Y)[test]
		Y_pred = clf.predict(np.array(X)[test])
		cm = confusion_matrix(Y_true, Y_pred)

		worksheet_write = book_write.add_sheet("confusion"+str(i))
		write_matrix(cm, worksheet_write)
		i += 1

def write_matrix(m, worksheet_write):
	for row in range(0,len(m)):
		for col in range(0, len(m[0])):
			write(row, col, m[row][col], worksheet_write)

def unbalanced_data_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]
	worksheet_write = book_write.add_sheet("weighted_classes")
	#AUTO adjust weights to be inversely proportional to class frequency
	clf = clf= svm.SVC(kernel='linear', class_weight='auto')

	row = 0
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	X = addFeatures(seqs, X, feature="ngram", n=2) #n=1-5
	val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
	write(row, 0, "ngrams, n=2", worksheet_write)
	write(row, 1, val, worksheet_write)

#Runs the experiment that scrapes uniprot to get the molecular function that we want.
def function_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]
	worksheet_write = book_write.add_sheet("functions")

	clf = clf= svm.SVC(kernel='linear', class_weight='auto')
	row = 0
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	'''scrape_info = scrape_list_ids(True)
	feature_dict = scrape_info[0]
	all_features = scrape_info[1]
	seq_ids = data[2]'''
	X = addFeatures(seqs, X, feature="ngram", n=3) 
	#X = addFeatures(seqs, X, feature="functions", ids=seq_ids, function_dict=feature_dict, all_functions=all_features) #n=1-5

	print 'done adding features'
	#print "all features: " + str(all_features) + "len: " + str(len(all_features))
	clf = clf= svm.SVC(kernel='linear', class_weight='auto')
	row = 0
	val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=5)
	write(row, 0, "function type", worksheet_write)
	write(row, 1, val, worksheet_write)
	print 'X: ' + str(X)
	print "VAL: " + str(val)

	print "END"


def write(x, y, value, worksheet):
	worksheet.write(x, y, value)

#poly_dimension_experiment()
#experiment_with_k()
#feature_experiment()
function_experiment()
#rbf_gridsearch_experiment()
#confusion_matrix_experiment()
#unbalanced_data_experiment()
book_write.save("svm_output.xls")


'''data = getData()
seqs = data[0]
X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
Y = data[1]

X = addFeatures(seqs, X, feature="ngram", n=2)
X = addFeatures(seqs, X, feature="mapseq", n=1)
clf = clf= svm.SVC(kernel='linear')
train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', [3], worksheet_write, 0)'''