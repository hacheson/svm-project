#!/usr/bin/env python
import sklearn
from sklearn import svm
from sklearn import datasets
def svm():
#	iris = datasets.load_iris()
#	digits = datasets.load_digits()

	from sklearn import svm
	X = [[0, 0], [1, 1]]
	y = [0, 1]
	clf = svm.SVC()
	clf.fit(X, y)  
	print clf
	print clf.predict([[2., 2.]])
	#SVC(C=1.0, cache_size=200, class_weight=None, coef0=0.0, degree=3,
	#gamma=0.0, kernel='rbf', max_iter=-1, probability=False, random_state=None,
	#shrinking=True, tol=0.001, verbose=False)




if __name__ == "__main__":
	svm()
 	#sys.exit(main())