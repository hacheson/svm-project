#!/usr/bin/python
import os
import sys
import utils

'''Hannah Acheson-Field, hacheson. 10-5-12. CS181.'''

def globalAlignment(seq1, seq2, matchScore, misMatchScore, gapScore):
	matchScore = int(matchScore)
	misMatchScore = int(misMatchScore)
	gapScore = int(gapScore)
	
	#File input
	file1 = open(fileA, 'r')
	file2 = open(fileB, 'r')
	readLines1 = file1.readlines()
	readLines2 = file2.readlines()
	linesA = utils.getArraysFromFile(readLines1)
	linesB = utils.getArraysFromFile(readLines2)
	
	'''Returns the optimal alignment of all sequences located in a fafsta file.
		ie if there are four sequences on one file and 5 on the other,
		this will find the alignments of all 20 combinations. This an array containing
		each of the sequences, the score, and the percent correct.'''
	ret = utils.initialize2DArray(len(linesA)*len(linesB), 5, 0)
	retIndex = 0
	idInfo1 = utils.getSequenceInfo(readLines1)
	idInfo2 = utils.getSequenceInfo(readLines2)
	for i in range(0, len(linesA)):
		for j in range(0, len(linesB)):
			info = getAlignment(linesA[i], linesB[j], matchScore, misMatchScore, gapScore)
			#alignments
			ret[retIndex][0] = info[0]
			ret[retIndex][1] = info[1]
			#score
			ret[retIndex][2] = info[2]
			#percent correct
			ret[retIndex][3] = idInfo1[i]
			ret[retIndex][4] = idInfo2[j]
			retIndex += 1
	file1.close()
	file2.close()
	return ret

'''Returns an optimal global alignment. This is where the bulk of the algorithm takes place.'''
def getAlignment(horzChars, vertChars, matchScore, misMatchScore, gapScore):
	scores = utils.initialize2DArray(len(horzChars)+1, len(vertChars)+1, 0)
	pointers = utils.initialize2DArray(len(horzChars)+1, len(vertChars)+1, [-float('inf'),-float('inf')])
	dirs = utils.initialize2DArray(len(horzChars)+1, len(vertChars)+1, -float('inf'))
	#print scores
	#sets the initial gap scores
	for i in range(0, len(horzChars)+1):
		scores[i][0] = i*gapScore
	for i in range(1, len(vertChars)+1):
		scores[0][i] = i*gapScore

	for i in range(1, len(horzChars)+1):
		for j in range(1, len(vertChars)+1):
			if horzChars[i-1] == vertChars[j-1]:
				scores[i][j] = utils.maxVal(scores[i-1][j-1]+matchScore, scores[i-1][j]+gapScore, scores[i][j-1]+gapScore)
				direc = utils.maxForPointer(scores[i-1][j-1]+matchScore, scores[i-1][j]+gapScore, scores[i][j-1]+gapScore)
			else:
				scores[i][j] = utils.maxVal(scores[i-1][j-1]+misMatchScore, scores[i-1][j]+gapScore, scores[i][j-1]+gapScore)
				direc = utils.maxForPointer(scores[i-1][j-1]+misMatchScore, scores[i-1][j]+gapScore, scores[i][j-1]+gapScore)
			
			if direc == 0:
				pointers[i][j] = [i-1,j-1]
				dirs[i][j] = direc
			elif direc == 1:
				pointers[i][j] = [i-1,j]
				dirs[i][j] = direc
			else:
				pointers[i][j] = [i, j-1]
				dirs[i][j] = direc

	maxScore = scores[len(horzChars)][len(vertChars)]
	pt = [len(horzChars),len(vertChars)]

	#traceback
	horzOut = []
	vertOut = []

	oldPt = pt
	while(pt[0]!= -float('inf') or pt[1]!= -float('inf')):
		if dirs[pt[0]][pt[1]] == 0:
			horzOut.insert(0, horzChars[pt[0]-1])
			vertOut.insert(0, vertChars[pt[1]-1])
		elif dirs[pt[0]][pt[1]] == 1:
			vertOut.insert(0, '-');
			horzOut.insert(0, horzChars[pt[0]-1])
		elif dirs[pt[0]][pt[1]] == 2:
			horzOut.insert(0, '-');
			vertOut.insert(0, vertChars[pt[1]-1])
		oldPt = pt
		pt = pointers[pt[0]][pt[1]]

	#printing indels on edge
	pt = oldPt
	if pt[0]!=0 or pt[1]!=0:
		#in vertical column
		if pt[1]!=0:
			numIndels = scores[0][pt[1]]/(gapScore)
			for i in range(0, numIndels):
				horzOut.insert(0, '-')
			while(pt[1]!=0):
				vertOut.insert(0, vertChars[pt[1]-1])
				pt[1] = pt[1]-1;
		elif pt[0] != 0:
			numIndels = scores[pt[0]][0]/(gapScore)
			for i in range(0, numIndels):
				vertOut.insert(0, '-');
			while(pt[0]!=0):
				horzOut.insert(0, horzChars[pt[0]-1])
				pt[0] = pt[0]-1;
			
	print "Score: " + str(maxScore)	
	#utils.printArray(horzOut)
	#utils.printArray(vertOut)
	#print '\n'
	return maxScore
	#return [''.join(horzOut), ''.join(vertOut), maxScore]	

if __name__ == '__main__':
	assert len(sys.argv) == 6
	globalAlignment(sys.argv[1],sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
