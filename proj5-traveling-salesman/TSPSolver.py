#!/usr/bin/python3

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

from copy import copy, deepcopy
import traceback
import time
import numpy as np
from TSPClasses import *
import heapq
import itertools
import traceback


class Stack(object):

    def __init__(self):
        self.items = []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def peek(self):
        return self.items[-1]

    def isEmpty(self):
        return len(self.items) == 0


class TSPSolver:
    def __init__(self, gui_view):
        self._scenario = None

    def setupWithScenario(self, scenario):
        self._scenario = scenario

    ''' <summary>
		This is the entry point for the default solver
		which just finds a valid random tour.  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of solution, 
		time spent to find solution, number of permutations tried during search, the 
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

    def defaultRandomTour(self, time_allowance=60.0):
        results = {}
        cities = self._scenario.getCities()
        ncities = len(cities)
        foundTour = False
        count = 0
        bssf = None
        start_time = time.time()
        while not foundTour and time.time() - start_time < time_allowance:
            # create a random permutation
            perm = np.random.permutation(ncities)
            route = []
            # Now build the route using the random permutation
            for i in range(ncities):
                route.append(cities[perm[i]])
            bssf = TSPSolution(route)
            count += 1
            if bssf.cost < np.inf:
                # Found a valid route
                foundTour = True

        end_time = time.time()
        results['cost'] = bssf.cost if foundTour else math.inf
        results['time'] = end_time - start_time
        results['count'] = count
        results['soln'] = bssf
        results['max'] = None
        results['total'] = None
        results['pruned'] = None
        return results

    ''' <summary>
		This is the entry point for the greedy solver, which you must implement for 
		the group project (but it is probably a good idea to just do it for the branch-and
		bound project as a way to get your feet wet).  Note this could be used to find your
		initial BSSF.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found, the best
		solution found, and three null values for fields not used for this 
		algorithm</returns> 
	'''

    def greedy(self, time_allowance=60.0):
        try:
            # initialization is O(n) space to hold all the cities and O(1) time
            bssf = None
            cities = self._scenario.getCities()
            self.num_cities = len(cities)
            solution_dict = {}
            start_time = time.time()
            # stay under the time limit
            while time.time() - start_time < time_allowance:
                # for each city as the starting node
                for index_loop in range(len(cities)):
                    print("Starting with city {}".format(index_loop))
                    city = cities[index_loop]
                    city_path = []
                    city_path.append(city)
                    to_visit = deepcopy(cities)
                    current_city = city
                    del to_visit[index_loop]
                    while len(to_visit):
                        city_costs = self.get_closest_cities(current_city, to_visit)
                        closest_city_tuple = city_costs[0]
                        closest_city_index = to_visit.index(closest_city_tuple[0])
                        closest_city = to_visit[closest_city_index]
                        if not self._scenario._edge_exists[current_city._index][closest_city._index]:
                            # should only happen if all costs to remaining cities are inf
                            break
                        del to_visit[closest_city_index]
                        city_path.append(closest_city)
                        current_city = closest_city
                    if len(to_visit):
                        continue
                    else:
                        bssf = TSPSolution(city_path)
                        end_time = time.time()
                        results = {}
                        results['cost'] = bssf.cost
                        results['time'] = end_time - start_time
                        results['count'] = None
                        results['soln'] = bssf
                        results['max'] = None
                        results['total'] = None
                        results['pruned'] = None
                        solution_dict[index_loop] = results
                        continue

                self.lowest_cost = float("inf")
                for key, solution in solution_dict.items():
                    print(key, solution["cost"])
                    if solution["cost"] < self.lowest_cost:
                        self.lowest_cost = solution["cost"]
                        lowest = solution

                print("BSSF is {}".format(lowest["soln"].cost))
                return lowest

        except Exception as e:
            print(e)
            traceback.print_exc()
            raise (e)

    def get_distance_matrix(self, city_list):
        ### create and initialize the matrix ###
        matrix = np.full((len(city_list), len(city_list)), fill_value=np.inf)
        for from_index, city in enumerate(city_list):
            for dest_index, dest_city in enumerate(city_list):
                if from_index == dest_index:
                    # already init to inf
                    continue
                dist = city.costTo(dest_city)
                matrix[from_index][dest_index] = dist

        return matrix

    ''' <summary>
		This is the entry point for the branch-and-bound algorithm that you will implement
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number solutions found during search (does
		not include the initial BSSF), the best solution found, and three more ints: 
		max queue size, total number of states created, and number of pruned states.</returns> 
	'''

    """
	ALGORITHM ANALYSIS:
	O(k*n^2) time where k is the number of states generated and checked
	O(k*n^2) space where k is the number of states generated and on the heap


	IMPLEMENTATION HELP:
	tuples are ordered as following in min heap
	0. value that will be used for sorting
	1. the current node to visit and expand
	2. the nodes that have yet to  be visited
	3. the current reduced matrix for that sub problem
	4. The order of the nodes visited
	5. The total cost of this subproblem 
	"""

    def branchAndBound(self, time_allowance=60.0):
        try:
            # initialization is O(n) space to hold all the cities and O(1) time
            bssf = self.greedy(time_allowance=time_allowance)['soln']
            # Random too slow
            # bssf2 = self.defaultRandomTour(time_allowance=time_allowance)['soln']
            # print("BSSF1 {} and BSSF2 {}".format(bssf.cost, bssf2.cost))
            # bssf = bssf if min(bssf.cost, bssf2.cost) == bssf.cost else bssf2
            print("Greedy and Random Algorithms found a BSSF of {}".format(bssf.cost))
            cities = self._scenario.getCities()
            self.cities = cities
            self.lowest_ave_cost = float("inf")
            self.num_cities = len(cities)
            heap = []
            bssf_updates = 0
            max_queue_len = 1
            total_states_created = 1
            pruned_subproblems = 0
            num_solutions = 0
            self.lowest_cost = bssf.cost
            # getting a recuced matrix is O(n^2) time and space
            initial_reduced_matrix, lower_bound = self.get_init_reduced_matrix(cities)
            print("Initial lower bound is {}".format(lower_bound))
            starting = tuple((lower_bound,
                              cities[0],
                              cities[1:],
                              initial_reduced_matrix,
                              [cities[0]._index],
                              lower_bound))
            heapq.heappush(heap, starting)
            time_start = time.time()
            # stay under the time limit
            while time.time() - time_start < time_allowance and len(heap):
                # find the best solutions
                # pop is O(logn) time to re-heapify and O(1) space
                next_to_try = heapq.heappop(heap)
                # check is O(1) time and space
                if next_to_try[5] < self.lowest_cost:
                    # for city in cities to visit still
                    # print("Exanding subproblem with cost {} and {} more cities to visit".format(next_to_try[5],
                    # 																			len(next_to_try[2])))
                    for city in next_to_try[2]:
                        # MAKE SURE PATH EXISTS
                        if self._scenario._edge_exists[next_to_try[1]._index][city._index]:
                            # getting reduced matrix is O(n^2) time and space
                            new_subproblem = self.get_reduced_matrix(city, next_to_try[3], next_to_try)
                            # if there are no more nodes to visit
                            if not len(new_subproblem[2]):
                                # this is O(n) to get the cities
                                route = self.index_to_cities(new_subproblem[4])
                                # checks and creation are O(1) time and space
                                bssf = TSPSolution(route)
                                if bssf.cost < self.lowest_cost:
                                    self.lowest_cost = min(bssf.cost, self.lowest_cost)
                                    bssf_updates += 1
                                num_solutions += 1
                                print("### Found a solution!!! Score: {} ####".format(self.lowest_cost))
                            # There are still nodes to visit
                            else:
                                # should we prune this node?
                                if new_subproblem[5] < self.lowest_cost:
                                    # It is O(logn) time to heapify after inserting
                                    # O(1) space
                                    heapq.heappush(heap, new_subproblem)
                                    total_states_created += 1
                                else:
                                    pruned_subproblems += 1
                                    total_states_created += 1


                else:
                    pruned_subproblems += 1
                    total_states_created += 1
                # checking is just comparing two numbers, O(1) time and space
                max_queue_len = max(len(heap), max_queue_len)

            # all of these are O(1) time and space to print and return
            end_time = time.time()
            print("BSSF updates: {}".format(bssf_updates))
            results = {}
            results['cost'] = self.lowest_cost
            results['time'] = end_time - time_start
            results['count'] = num_solutions
            results['soln'] = bssf
            results['max'] = max_queue_len
            print("The max length of the queue is {}".format(max_queue_len))
            results['total'] = total_states_created
            print("The total number of nodes created was {}".format(total_states_created))
            results['pruned'] = pruned_subproblems
            print("The total number of nodes pruned was {}".format(pruned_subproblems))

            return results
        except Exception as e:
            # error catching so I can see errors
            print(e)
            print(traceback.format_exc())
            raise (e)

    """
	This involves only a lookup and a divide so it is O(1) time and space to get
	the length of the array and divide
	"""

    def get_value(self, score, cities_visited):
        if len(cities_visited):
            # note I tried a lot of these to get varying results
            return score / len(cities_visited)
        else:
            # we will add this to bssf anyways
            return 0

    """
	This function is O(n^2) time and space to visit and update every cell and return an array of 
	size O(n^2). More details are in the comments.
	"""

    def get_reduced_matrix(self, next_city_to_visit, matrix, given_tuple):
        # copy the tuple so it doesn't make changes to the main one
        # these are O(1) time and O(n^2) space to copy the matrix
        tuple_of_items = deepcopy(given_tuple)
        matrix = matrix.copy()
        sum_to_reduce = 0
        initial_cost_to_city = matrix[tuple_of_items[1]._index][next_city_to_visit._index]

        # since we're visiting the row and column, make them inf
        # O(1) time and space
        matrix[tuple_of_items[1]._index] = np.inf
        matrix[:, next_city_to_visit._index] = np.inf
        matrix[next_city_to_visit._index][tuple_of_items[1]._index] = np.inf

        # The following code runs over all the rows and then all the columns.  This is
        # O(n^2) where n is the number of cities since it updates every part of the matrix
        # it is O(1) space since it does it with the same matrix

        # do rows and get mins
        for row in range(matrix.shape[0]):
            min_for_row = np.min(matrix[row])
            if np.isinf(min_for_row):
                continue
            matrix[row] = matrix[row] - min_for_row
            sum_to_reduce += min_for_row

        # get reduced columns now
        for col in range(matrix.shape[1]):
            min_for_col = np.min(matrix[:, col])
            if np.isinf(min_for_col):
                continue
            matrix[:, col] = matrix[:, col] - min_for_col
            sum_to_reduce += min_for_col

        # This is O(n) time to visit all the cities and O(1) space to delete it
        cities_to_still_visit = tuple_of_items[2]
        cities_to_still_visit = self.delete_city_by_id(cities_to_still_visit, next_city_to_visit)

        new_cost = tuple_of_items[5] + initial_cost_to_city + sum_to_reduce

        ### Create the tuple to store in the queue: this is O(n^2) space to hold the matrix
        return_tuple = ((self.get_value(tuple_of_items[5], tuple_of_items[4]),  # the value for queue
                         next_city_to_visit,  # the city object that was visited now
                         cities_to_still_visit,  # the nodes to still visit
                         matrix,  # the new distance matrix
                         tuple_of_items[4] + [next_city_to_visit._index],  # the order of nodes visited
                         new_cost))  # the new cost

        return return_tuple

    """
	This is a O(n^2) algorithm to get the distance from all nodes to every other node
	Also O(n^2) space to hold the matrix
	"""

    def get_init_reduced_matrix(self, city_list):
        ### create and initialize the matrix ###
        matrix = np.full((len(city_list), len(city_list)), fill_value=np.inf)
        for from_index, city in enumerate(city_list):
            for dest_index, dest_city in enumerate(city_list):
                if from_index == dest_index:
                    # already init to inf
                    continue
                dist = city.costTo(dest_city)
                matrix[from_index][dest_index] = dist

        # The following options are O(n^2) time to operate on every cell
        # it is O(1) space
        sum_to_reduce = 0
        # do rows and get mins
        for row in range(matrix.shape[0]):
            min_for_row = np.min(matrix[row])
            matrix[row] = matrix[row] - min_for_row
            sum_to_reduce += min_for_row

        # get reduced columns now
        for col in range(matrix.shape[1]):
            min_for_col = np.min(matrix[:, col])
            matrix[:, col] = matrix[:, col] - min_for_col
            sum_to_reduce += min_for_col

        return matrix, sum_to_reduce

    """
	This is O(1) time and space to do a multiplication and a subtraction
	This function is seperate so that it can be changed easily
	"""

    def get_cost_from_value(self, new_subproblem):
        value = new_subproblem[0]
        num_cities_left = len(new_subproblem[2])
        return value  ## / np.square(self.num_cities - num_cities_left)

    """"
	This function turns all the indexes into their City class counterpart
	It is O(n) time and space to iterate through and hold all the cities
	"""

    def index_to_cities(self, city_indices):
        city_list = []
        for city_int in city_indices:
            city_list.append(self.cities[city_int])

        return city_list

    """
	This is just a print debug function.  O(1) time and space
	"""

    def print_subprob(self, subprob):
        print("Node: {}, current_cost: {}".format(subprob[1]._index, subprob[5]))

    """
	This function goes through the list and delete the node based on it's ID
	It is O(n) time to go through every node worst case and O(1) space
	"""

    def delete_city_by_id(self, cities_to_still_visit, next_city_to_visit):
        for index, city in enumerate(cities_to_still_visit):
            if city._index == next_city_to_visit._index:
                delete_index = index
                break

        del cities_to_still_visit[delete_index]
        return cities_to_still_visit

    """
	Used for greedy algorithm
	This is O(n) time and space to find the closest city and to store the list and find the min
	"""

    def get_closest_cities(self, city, city_list):
        cost = {}
        for city_to_visit in city_list:
            cost[city_to_visit] = city.costTo(city_to_visit)

        sorted_x = sorted(cost.items(), key=lambda kv: kv[1])
        # print("Closest length is {}".format(sorted_x[0][1]))
        return sorted_x

    ''' <summary>
		This is the entry point for the algorithm you'll write for your group project.
		</summary>
		<returns>results dictionary for GUI that contains three ints: cost of best solution, 
		time spent to find best solution, total number of solutions found during search, the 
		best solution found.  You may use the other three field however you like.
		algorithm</returns> 
	'''
    from collections import Counter
    def fancy(self, time_allowance=60.0):
        try:
            start_node_num = 0
            cost = {}
            city_list = self._scenario.getCities()
            start_city = city_list[start_node_num]
            start_time = time.time()
            for city_to_visit in city_list:
                cost[frozenset([city_to_visit._index])] = {"cost": start_city.costTo(city_to_visit),
                                                           "prev": city_to_visit._index}

            cities_without_start = [number for number in range(len(city_list)) if number != start_node_num]
            while time.time() - start_time < time_allowance:
                for subset_length in range(2, len(city_list)):
                    combinations_list = list(itertools.combinations(cities_without_start, subset_length))
                    for combination in combinations_list:
                        min_value = float("inf")
                        min_prev = None
                        for city_index, city_num in enumerate(combination):
                            subset = list(combination)
                            del subset[city_index]
                            lookup = cost[frozenset(subset)]
                            value = lookup["cost"] + city_list[lookup["prev"]].costTo(city_list[city_num])
                            if value <= min_value:
                                min_value = value
                                min_prev = city_num
                        cost[frozenset(combination)] = {"cost": min_value, "prev": min_prev}

                # do last iteration
                min_value = float("inf")
                min_prev = None
                # try to see distance to initial city from all other cities
                for city_index, city_num in enumerate(cities_without_start):
                    subset = list(cities_without_start)
                    del subset[city_index]
                    lookup = cost[frozenset(subset)]
                    # final cost is cost from city_num as 2nd to last coming from the rest of the subset
                    # and going to the zero-th node
                    value = lookup["cost"] + city_list[lookup["prev"]].costTo(city_list[city_num]) + \
                            city_list[city_num].costTo(city_list[0])
                    if value <= min_value:
                        min_value = value
                        min_prev = city_num
                cost[frozenset(city_list)] = {"cost": min_value, "prev": min_prev}

                # trace pointer back
                route = []
                # this is the last node
                current_prev_pointer = cost[frozenset(city_list)]["prev"]
                current_cities_set = set(cities_without_start)
                # go until 2nd to last node (otherwise this method doesn't work)
                while len(current_cities_set) > 1:
                    # insert node into route list
                    route.insert(0, city_list[current_prev_pointer])
                    # get set without that city
                    current_cities_set = current_cities_set - set([current_prev_pointer])
                    current_prev_pointer = cost[frozenset(current_cities_set)]["prev"]

                # add last element of set - 2nd to beginning node
                route.insert(0, city_list[current_cities_set.pop()])
                # add start city
                route.insert(0, city_list[start_node_num])
                # return solution
                bssf = TSPSolution(route)
                end_time = time.time()
                results = {}
                results['cost'] = bssf.cost
                results['time'] = end_time - start_time
                results['count'] = None
                results['soln'] = bssf
                results['max'] = None
                results['total'] = None
                results['pruned'] = None
                return results

        except Exception as e:
            # error catching so I can see errors
            print(e)
            print(traceback.format_exc())
            raise (e)
