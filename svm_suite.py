from svm import *


#example
#opens workbooks
book = open_workbook('Adomain_Substrate.xls')
worksheet = book.sheet_by_name('Adomain_Substrate')
book_write = Workbook()
worksheet_write = book_write.add_sheet("Kernels")

def poly_dimension_experiment():
	data = getData()
	seqs = data[0]
	X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
	Y = data[1]

	X = addFeatures(seqs, X, feature="ngram", n=2)
	X = addFeatures(seqs, X, feature="mapseq", n=1)
	for i in range(2, 10):
		clf = clf= svm.SVC(kernel='poly', degree=i)
		val = train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', [3])
		write(0, i, val, worksheet_write)	
def write(x, y, value, worksheet):
	worksheet_write.write(x, y, value)

poly_dimension_experiment()
book_write.save("svm_output.xls")

'''data = getData()
seqs = data[0]
X = [ [] for _ in range(0, len(seqs))] # Empty feature vector for every sequence
Y = data[1]

X = addFeatures(seqs, X, feature="ngram", n=2)
X = addFeatures(seqs, X, feature="mapseq", n=1)
clf = clf= svm.SVC(kernel='linear')
train_test_SVM(np.array(X), np.array(Y), clf, 'k_fold', [3], worksheet_write, 0)'''