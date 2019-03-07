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
		self.matrix = np.zeros(len(sequences[0]), len(sequences[1]))
		results = []
		for i in range(len(sequences)):
			jresults = []
			for j in range(len(sequences)):

				if(j < i):
					s = {}
				else:
					score = self.get_gene_sequences(i, j)
					alignment1 = 'abc-easy  DEBUG:(seq{}, {} chars,align_len={}{})'.format(i+1,
						len(sequences[i]), align_length, ',BANDED' if banded else '')
					alignment2 = 'as-123--  DEBUG:(seq{}, {} chars,align_len={}{})'.format(j+1,
						len(sequences[j]), align_length, ',BANDED' if banded else '')
					s = {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
					table.item(i,j).setText('{}'.format(int(score) if score != math.inf else score))
					table.repaint()	
				jresults.append(s)
			results.append(jresults)
		return results

	def get_gene_sequences(self, i: int, j: int):
		top = float("inf")
		diag = float("inf")
		left = float("inf")

		# Get previous info
		if j - 1 > 0:
			top = self.matrix[i][j-1]
		if i - 1 > 0:
			left = self.matrix[i-1][j]
		if i - 1 > 0 and j - 1 > 0:
			diag = self.matrix[i - 1][j - 1]

		diag_value = SUB if self.string1[i] != self.string2[j] else MATCH
		choices = {top + INDEL: "top", left + INDEL: "left", diag + diag_value: "diag"}
		min_score = min(list(choices.keys()))
		# return min_score, choices[min_score]
		return 1 + choices




