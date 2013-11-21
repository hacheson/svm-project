from svm import *
from xlrd import open_workbook
from xlwt import *


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
	col = 0

	#AAcounts
	X = addFeatures(seqs, X, feature="AAcounts")
	val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
	write(col, 1, val, worksheet_write)
	write(col, 0, "AAcounts", worksheet_write)
	col+=1

	#AAncounts
	for i in range(1,6):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="AAncounts", n=i) #n=1-5
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(col, 0, "AAncounts, n=" + str(i), worksheet_write)
		write(col, 1, val, worksheet_write)
		col+=1
	#mapseq
	for i in range(1,6):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="mapseq", n=i) #n=1-5
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(col, 0, "mapseq, n=" + str(i), worksheet_write)
		write(col, 1, val, worksheet_write)
		col+=1

	#ngrams
	for i in range(1,4):
		X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
		X = addFeatures(seqs, X, feature="ngram", n=i) #n=1-5
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'strat_k_fold', k=4)
		write(col, 0, "ngrams, n=" + str(i), worksheet_write)
		write(col, 1, val, worksheet_write)
		col+=1

def write(x, y, value, worksheet):
	worksheet.write(x, y, value)

#poly_dimension_experiment()
#experiment_with_k()
feature_experiment()
book_write.save("svm_output.xls")


'''data = getData()
seqs = data[0]
X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
Y = data[1]

X = addFeatures(seqs, X, feature="ngram", n=2)
X = addFeatures(seqs, X, feature="mapseq", n=1)
clf = clf= svm.SVC(kernel='linear')
train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', [3], worksheet_write, 0)'''