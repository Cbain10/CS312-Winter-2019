#!/usr/bin/python3
# why the shebang here, when it's imported?  Can't really be used stand alone, right?  And fermat.py didn't have one...
# this is 4-5 seconds slower on 1000000 points than Ryan's desktop...  Why?
import time
import math
import itertools
from functools import partial

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF, QThread, pyqtSignal
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF, QThread, pyqtSignal
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))


"""
Since the QPointF class does not define the less than operator, 
I have defined it to sort via the x-axis.  This makes the sorting operation O(nlogn) since
that's what python uses (it uses Timsort) and the comparison is O(1): 
https://stackoverflow.com/questions/14434490/what-is-the-complexity-of-this-python-sort-method
https://en.wikipedia.org/wiki/Timsort
Space complexity is O(n)
"""
def less_than(self, other):
	return self.x() < other.x()
QPointF.__lt__ = less_than

"""
for closewiseangle_and_distance
Time complexity: O(1) since it is just mathematical operators
Space complexity: O(1) no extra space
"""
def clockwiseangle_and_distance(point, center_x, center_y, anchor="ne"):
	return math.atan2(point.x() - center_x, point.y() - center_y)

""""
this function is O(n) time complexity since it goes through every point once
and the divides and references are constant time
Space complexity is O(n) since it stores the points once again and makes the copy of n points
"""
def get_center(points_list):
	x = [p.x() for p in points_list]
	y = [p.y() for p in points_list]
	return sum(x) / len(points_list), sum(y) / len(points_list)

"""
This function is O(n) since it goes through every point and compares them with the max
Space complexity is O(1)
"""
def get_closest_x(points, side_to_find="left"):
	closest_index = 0
	if side_to_find == "left":
		for index, point in enumerate(points):
			if point > points[closest_index]:
				closest_index = index
	else:
		for index, point in enumerate(points):
			if point < points[closest_index]:
				closest_index = index
	return closest_index

class ConvexHullSolverThread(QThread):
	def __init__( self, unsorted_points,demo):
		self.points = unsorted_points					
		self.pause = demo
		QThread.__init__(self)

	def __del__(self):
		self.wait()

	show_hull = pyqtSignal(list,tuple)
	display_text = pyqtSignal(str)

# some additional thread signals you can implement and use for debugging, if you like
	show_tangent = pyqtSignal(list,tuple)
	erase_hull = pyqtSignal(list)
	erase_tangent = pyqtSignal(list)

	"""
	Divide and conquer is T(n) = 2T(n/2) + O(merge = n)
	Using the masters theorem this algorithm is O(nlogn)
	space complexity is log(n) stack calls of size n/2 so it is O(nlogn)
	"""
	def divide_and_conquer(self, points):
		# base case -> now go up the stack and merge
		if len(points) == 1:
			# this is O(1)
			return points

		# not base case, divide more
		# each call of T(n) = 2T(n/2) + O(merge)
		left_half = self.divide_and_conquer(points[0: len(points) // 2])
		right_half = self.divide_and_conquer(points[len(points) // 2:])

		# now that we have 2 hulls we can merge them
		# this step is O(n) time and O(n) space
		return self.merge(left_half, right_half)

	# return slope of two points
	"""
	get_slope is is O(1) since it just two subtractions and a divide
	Space complexity is O(1)
	"""
	@staticmethod
	def get_slope(left_point, right_point):
		return (right_point.y() - left_point.y()) / (right_point.x() - left_point.x())

	"""
	Time complexity is O(n) since we go through a variety of O(n) steps to find the tangents.
	Space complexity is O(n) since it creates a new array with "sorted" and get_center is O(n)
	"""
	def merge(self, left_half, right_half):
		if len(left_half) + len(right_half) < 4:
			# 3 or less points, return
			# this function is O(n) described above function
			center_x, center_y = get_center(left_half + right_half)
			# sorting is O(nlogn) as explained above but since this is the base case it is constant
			# #since at max it is 3 points
			full = sorted(left_half + right_half, key=partial(clockwiseangle_and_distance, center_x = center_x,
															  center_y = center_y))
			return full

		# get the index of the points closest to each other on the x-axis
		# these steps are O(n)
		left_start = get_closest_x(left_half, "left")
		right_start = get_closest_x(right_half, "right")

		# these steps are o(n) time worst case
		left_upper_tangent, right_upper_tangent = self.get_upper_tangent(left_start, left_half, right_start, right_half)
		# lower tangent is the reverse of upper tangent so switch right and left
		right_lower_tangent, left_lower_tangent = self.get_upper_tangent(right_start, right_half, left_start, left_half)
		# remove unneeded points - by creating a circular list
		# going through the points is at most O(n)
		full_list = []
		full_list.append(left_half[left_upper_tangent])
		index = right_upper_tangent
		while right_half[index % len(right_half)] != right_half[right_lower_tangent]:
			full_list.append(right_half[index % len(right_half)])
			index += 1
		full_list.append(right_half[right_lower_tangent])
		index = left_lower_tangent
		while left_half[index % len(left_half)] != left_half[left_upper_tangent]:
			full_list.append(left_half[index % len(left_half)])
			index += 1

		self.plot_points(full_list)

		return full_list

	"""
	get_upper_tangent: worst case we go through every node in the while loop which is O(n) worst case
	Space complexity is O(1)
	"""
	def get_upper_tangent(self, left_start, left_half, right_start, right_half):
		left_index = left_start
		right_index = right_start
		previous_tangent = tuple([left_index, right_index])
		while True:
			slope = self.get_slope(left_half[left_index], right_half[right_index])
			# try to move the slope again
			left_slope_decreasing = True
			# move left counter-clockwise as long as it makes left turn on LEFT HULL
			while left_slope_decreasing:
				next_iteration_slope = self.get_slope(left_half[(left_index - 1) % len(left_half)],
													  right_half[right_index])
				if next_iteration_slope < slope:
					slope = next_iteration_slope
					left_index -= 1
					# this step can only happen once
					if left_index < 0:
						left_index = len(left_half) - 1
					left_slope_decreasing = True
				else:
					left_slope_decreasing = False

			# try to move the slope again
			right_slope_increasing = True
			# move right clockwise as long as it makes right turn on RIGHT HULL
			while right_slope_increasing:
				next_iteration_slope = self.get_slope(left_half[left_index],
													  right_half[(right_index + 1) % len(right_half)])
				if next_iteration_slope > slope:
					slope = next_iteration_slope
					right_index += 1
					# this step can only happen once
					if right_index >= len(right_half):
						right_index %= len(right_half)
					right_slope_increasing = True
				else:
					right_slope_increasing = False


			current_tangent = tuple([left_index, right_index])
			# if the loop has run once and hasn't changed, quit
			if previous_tangent == current_tangent:
				break

			previous_tangent = current_tangent

		#print("Final Tangent points are {} and {}".format(left_half[left_index], right_half[right_index]))

		return left_index, right_index

	"""
	The run itself is the same as divide_and_conquer so it is 
	Time complexity: O(nlogn)
	Space complexity: O(nlogn)
	"""
	def run( self):
		assert( type(self.points) == list and type(self.points[0]) == QPointF )

		n = len(self.points)
		print('Computing Hull for set of {} points'.format(n) )

		t1 = time.time()
		# This step is O(nlogn)
		self.points.sort()
		t2 = time.time()
		print('Time Elapsed (Sorting): {:3.3f} sec'.format(t2-t1))

		t3 = time.time()
		self.lines = []
		# this step is O(nlogn)
		hull_pts = self.divide_and_conquer(self.points)
		t4 = time.time()

		# this is O(n) max to get all the points
		polygon = [QLineF(hull_pts[i], hull_pts[(i + 1) % len(hull_pts)]) for i in range(len(hull_pts))]
		# when passing lines to the display, pass a list of QLineF objects.  Each QLineF
		# object can be created with two QPointF objects corresponding to the endpoints
		assert (type(polygon) == list and type(polygon[0]) == QLineF)
		# send a signal to the GUI thread with the hull and its color
		self.show_hull.emit(polygon, (255, 0, 0))
		# TODO: PASS THE CONVEX HULL LINES BACK TO THE GUI FOR DISPLAY

		# send a signal to the GUI thread with the time used to compute the hull
		self.display_text.emit('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4-t3))
		print('Time Elapsed (Convex Hull): {:3.3f} sec'.format(t4-t3))
			

	# debug function, O(n)
	def plot_points(self, hull_pts):
		polygon = [QLineF(hull_pts[i], hull_pts[(i + 1) % len(hull_pts)]) for i in range(len(hull_pts))]
		# when passing lines to the display, pass a list of QLineF objects.  Each QLineF
		# object can be created with two QPointF objects corresponding to the endpoints
		assert (type(polygon) == list and type(polygon[0]) == QLineF)
		# send a signal to the GUI thread with the hull and its color
		self.show_hull.emit(polygon, (255, 0, 0))

