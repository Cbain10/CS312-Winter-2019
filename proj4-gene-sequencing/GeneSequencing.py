#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import numpy as np

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class Cell:

	def __init__(self, cost=0, prev=None):
		self.cost = cost
		self.prev = prev


class GeneSequencing:

	def __init__( self ):
		self.matrix = np.ndarray

	
# This is the method called by the GUI.  _sequences_ is a list of the ten sequences, _table_ is a
# handle to the GUI so it can be updated as you find results, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you 
# how many base pairs to use in computing the alignment

	def align( self, sequences, table, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

		results = []
		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):

				if(j < i):
					s = {}
				else:
					try:
						print("On {}, {}".format(i, j))
						score, alignment1, alignment2 = self.get_gene_sequences(sequences[i], sequences[j],
																				align_length, banded)
					except Exception as e:
						# to debug
						print(e)
						quit(1)
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()	
				jresults.append(s)
			results.append(jresults)
		return results


	def get_gene_sequences(self, string1, string2, align_length, banded=False, k = 7):
		if string1 == string2:
			return 0, string1, string2
		else:
			string1 = string1[:align_length]
			string2 = string2[:align_length]

			if not banded:
				# init the matrices.  This takes O(M*N) twice which is still O(M*N) time and space
				self.matrix = np.full(shape=(len(string1)+ 1, len(string2) + 1), fill_value=float("inf"))
				self.matrix[0][0] = 0
				self.prev_matrix = np.empty(shape=(len(string1) +1, len(string2)+1), dtype=str)
			else:
				# init the matrices.  This takes k * O(N) twice which is still O(k*N) time and space
				self.matrix = np.full(shape=(len(string1) + 1, k), fill_value=float("inf"))
				self.matrix[0][0] = 0
				self.prev_matrix = np.empty(shape=(len(string1) + 1, k), dtype=str)

		# Outer loop runs M times
		shift_threshold = (k // 2)
		for i in range( len(string1)+1):
			# Inner loop runs N times if not-banded
			for j in range(len(string2) + 1):
				if not banded and i > j:
					continue
				# if banded run 7 * M times max, or O(M)
				if banded and abs(j - i) > (k // 2):
					continue
				top = float("inf")
				diag = float("inf")
				left = float("inf")
				new_option_flag = False
				# Get previous info
				if not banded:
					if i - 1 >= 0:
						new_option_flag = True
						top = self.matrix[i - 1][j]
					if j - 1 >= 0:
						new_option_flag = True
						left = self.matrix[i][j - 1]
					if i - 1 >= 0 and j - 1 >= 0:
						new_option_flag = True
						diag = self.matrix[i - 1][j - 1]
				else:
					if i > shift_threshold:
						shift = i - shift_threshold
					else:
						shift = 0
					if i - 1 >= 0 and j - shift + 1 < k:
						new_option_flag = True
						top = self.matrix[i - 1][j + 1 - shift]
					if j - shift - 1 >= 0:
						new_option_flag = True
						left = self.matrix[i][j - 1 - shift]
					if i - 1 >= 0:
						new_option_flag = True
						diag = self.matrix[i - 1][j - shift]
				print(i, j, "After")
				if new_option_flag:
					# edge cases, literally
					if not i or (not j and not banded):
						diag_value = float('inf')
					else:
						diag_value = SUB if string1[i - 1] != string2[j - 1] else MATCH
					# Since this compares three values to get min it is O(1) for this and other operations
					choices = {top + INDEL: "top", left + INDEL: "left", diag + diag_value: "diag"}
					min_score = min(list(choices.keys()))
					if not banded:
						self.matrix[i][j] = min_score
						self.prev_matrix[i][j] = choices[min_score]
					else:
						self.matrix[i][j - shift] = min_score
						self.prev_matrix[i][j - shift] = choices[min_score]

		if not banded:
			optimal_value = self.matrix[len(string1),  len(string2)]
			alignment1, alignment2 = self.get_alignments(string1, string2)
		else:
			optimal_value = self.matrix[len(string1),  shift_threshold + 1]
			alignment1, alignment2 = self.get_alignments_banded(string1, string2, k)

		# This step is the max of O(M) and O(N)
		return optimal_value, alignment1, alignment2

	"""
	Time complexity: goes backwards to generate the string.  Max of O(N) or O(M).
	Space complexity: max of O(M) or O(N) to hold a string
	"""
	def get_alignments(self, string1, string2):
		i = len(string1)
		j = len(string2)
		alignment1 = string1
		alignment2 = string2
		new = ""
		while True:
			prev = self.prev_matrix[i][j]
			if prev == "":
				break
			if prev == "t":
				new = "-" + new
				i = i - 1
				pass
			elif prev == "l":
				new = "-" + new
				j = j - 1
				pass
			elif prev == "d":
				new = string1[i - 1] + new
				i = i - 1
				j = j - 1
				pass
			else:
				print('Incorrect Previous Type. Bug!')
				raise Exception
		return alignment1, new

	"""
		Time complexity: goes backwards to generate the string.  Max of O(N) or O(M).
		Space complexity: max of O(M) or O(N) to hold a string
	"""
	def get_alignments_banded(self, string1, string2, k):
		i = len(string1)
		j = (k // 2) + 1

		new1 = ""
		new2 = ""
		while True:
			prev = self.prev_matrix[i][j]
			if prev == "":
				break
			if prev == "t":
				new2 = "-" + new2
				i = i - 1
				j = j + 1
				pass
			elif prev == "l":
				new1 = "-" + new1
				j = j - 1
				pass
			elif prev == "d":
				new1 = string1[i - 1] + new1
				new2 = string2[j - 1] + new2
				i = i - 1
				pass
			else:
				print('Incorrect Previous Type. Bug!')
				raise Exception
		return new1, new2







