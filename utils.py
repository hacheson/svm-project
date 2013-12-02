'''Utils functions for sequence alignments. 
	Hannah Acheson-Field, hacheson. 10-5-12'''


'''Returns the max of three values. Used instead of built in functionality of python
	to be sure which value is being returned.'''
def maxVal(a, b, c):
	if a >= b and a >= c:
		return a
	elif b>=a and b>=c:
		return b
	else:
		return c

'''Returns the max of two values'''
def maxValTwo(a, b):
	if a>= b:
		return a
	else:
		return b

'''Returns a coding depending on which direction: 0=diagonal, 1=left, 2=up'''
def maxForPointer(a, b, c):
	if a >= b and a >= c:
		return 0
	elif b>=a and b>=c:
		return 1
	else:
		return 2

'''Initializes 2D arrays: dimensions x and y, and initial value.
	Initialize to x+1 or y+1 becuase of initialization'''
def initialize2DArray(x, y, val):
	ret = []
	for i in range(0, x):
		ret.append([])
		for j in range(0, y):
			ret[i].append(val)
	return ret	

'''Prints out the arrays'''
def printArray(array):
	out = ""
	for i in range(0, len(array)):
		out += array[i]
	print out

'''Returns a 2D array of sequences broken down by lines starting with a >'''
def getArraysFromFile(lines):
	sequences = []
	for i in range(0, len(lines)):
		if lines[i].startswith('>'):
			sequences.append(createArrayAfterLine(i+1, lines))
	return sequences

'''Does similar process as getArraysFromFile, but this returns an array of the info of the strings'''	
def getSequenceInfo(lines):
	sequences = []
	info = []
	for i in range(0, len(lines)):
		if lines[i].startswith('>'):
			info.append(lines[i])
	return info

'''Creates an array of strings after a certain line stopping if a line begins with a >'''				
def createArrayAfterLine(index, lines):
	ret = []
	while lines[index].startswith('>') == False:
		ret.append(lines[index])
		index += 1
		if not index<len(lines):
			break
	return getChars(ret)

'''Returns an array of characters given an array of strings'''
def getChars(array):
	ret = []
	for i in range(0, len(array)):
		string = array[i].rstrip('\n')
		chars = list(string)
		for j in range(0, len(chars)):
			ret.append(chars[j])
	return ret

'''Returns the max pointer for four values.'''
def maxForPointerFour(a, b, c, d):
	if a >= b and a >= c and a >= d:
		return 0
	elif b>=a and b>=c and b>=d:
		return 1
	elif c>=a and c>=b and c>=d:
		return 2
	else:
		return 3

'''Returns the max of four values.'''
def maxValFour(a, b, c, d):
	if a >= b and a >= c and a >= d:
		return a
	elif b>=a and b>=c and b>=d:
		return b
	elif c>=a and c>=b and c>=d:
		return c
	else:
		return d
