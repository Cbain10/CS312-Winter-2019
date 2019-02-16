#!/usr/bin/python3
import traceback

from CS312Graph import *
import time


class NetworkRoutingSolver:
    def __init__( self ):
        pass

    def initializeNetwork( self, network ):
        assert( type(network) == CS312Graph )
        self.network = network

    """
    worst case we go through all V nodes to get there.  For each node worst case we check each edge attatched.
    Thus, this worst case O(V + E) time and O(V + E) space to hold all those nodes and edges.
    """
    def getShortestPath( self, destIndex ):
        self.dest = destIndex
        # init O(1) time and space
        path_edges = []
        total_length = 0
        nodes = self.network.nodes
        prev = self.node_info["prev"]
        cost = self.node_info["cost"]
        current_node_index = destIndex
        # worst case we go through all V nodes to get there.  For each node worst case we check each edge attatched.
        # thus this worst case O(V + E) time and O(V + E) space to hold all those nodes and edges.
        while current_node_index != self.source:
            prev_node = prev[current_node_index]
            if prev_node is None:
                print("There is no previous node.  Quitting")
                return {'cost': float("inf"), 'path': []}
            # check all the edges worst case
            edge_length = self.get_node_length(prev_node, current_node_index)
            if edge_length is None:
                print("There is no edge.  Quitting")
                return {'cost': float("inf"), 'path': []}
            path_edges.append((nodes[current_node_index].loc, nodes[prev_node].loc, '{:.0f}'.format(edge_length)))
            total_length += edge_length
            current_node_index = prev_node
        assert (int(total_length) == int(cost[destIndex]))
        return {'cost': cost[destIndex], 'path': path_edges}

    """
    This function is the same as dijkstras.  Since it is such a long explanation, see dijkstras for more info
    """
    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        self.node_info = self.dijkstras(self.network.nodes, self.source, use_heap)
        t2 = time.time()
        return (t2-t1)

    """
    This checks at most all the edges for a node which is O(1) = 3 worst case. It is O(1) space since it 
    doesn't store anything
    """
    def get_node_length(self, current_node_index, prev_node):
        for edge in self.network.nodes[current_node_index].neighbors:
            if edge.dest.node_id == prev_node:
                return edge.length

    """
    Depending on the implementation (list or heap) we see a time complexity of O(V^2 + E) or O((V+E)logV) 
    time respectively.
    Depending on the implementation (list or heap) we see a space complexity of O(V) space for both.
    
    See details in function or in paper report.
    """
    def dijkstras(self, graph, start_node, use_heap):
        try:
            # initialize -> O(1)
            dist = {}
            prev = {}
            length_dict = {}
            # for loop is O(V) time and space
            for index, node in enumerate(graph):
                dist[node.node_id] = float('inf')
                prev[node.node_id] = None

            dist[start_node] = 0
            # set up data structures
            if use_heap:
                priority = MinHeap(dist)  # O(VlogV) time to insert V nodes at log(V) cost.
            else:
                priority = Queue(start_node, dist.keys())  # O(V) time to insert V nodes at O(1) cost.

            iteration = 0
            # Checking size is O(1) worst case for both implementations
            while priority.size() != 0:
                # delete min is O(logV) for heap and O(V) for list.  Space is O(1) for both
                # thus sine this is called V times it is O(VlogV) for heap and O(V^2) for list
                u = priority.delete_min()
                # For every edge: all e in E
                for edge in graph[u].neighbors:
                    if edge.length + dist[u] < dist[edge.dest.node_id]:
                        dist[edge.dest.node_id] = dist[u] + edge.length
                        prev[edge.dest.node_id] = u
                        # decrease key is called E times worst case
                        # for heap it is O(logV) and for list O(1).
                        # thus it is O(ElogV) for the heap and O(E) for the list
                        priority.decrease_key(edge.dest.node_id, dist[edge.dest.node_id])
                iteration += 1
            return {"cost": dist, "prev": prev}

        except Exception as e:
            print(e)
            traceback.print_tb(e.__traceback__)


class Queue:

    """
    creating the full list is n times of inserting, which is O(V)
    Space is O(V) since it stores all the nodes
    """
    def __init__(self, start_node, list_of_all_nodes):
        # self.queue = [float("inf")] * len(list_of_all_nodes)
        self.queue = []
        # insert for loop -> n times O(1) or O(n)
        for node_id in list_of_all_nodes:
            self.insert(node_id, float("inf"))

        # for this lab set the start_node to zero distance
        self.queue[start_node] = 0
        self.length = len(self.queue)

    """
    O(1) time and space
    Note: this is never used in this data structure -> only outside
    """
    def size(self):
        return self.length

    """
    It is O(V) time since it has to find the min value and return it.
    It is O(1) Space since it doesn't store anything
    """
    def delete_min(self):
        # this is O(n) since it takes n to find the min and n to find the index
        min_index = self.queue.index(min(x for x in self.queue if x is not None))
        self.queue[min_index] = None
        self.length -= 1
        return int(min_index)

    """
    It is O(1) time since it just changes the value of an index
    It is O(1) Space since it doesn't store anything
    """
    def decrease_key(self, key, value):
        # changing the value is O(1) for lookups
        self.queue[key] = value
        return

    """
    It is O(1) time since it just appends to the end of a list and changes it's value
    It is O(1) Space since it doesn't store anything other than the one value
    """
    def insert(self, key, value):
        self.queue.append(key)
        self.queue[key] = value
        self.length += 1
        return

class MinHeap:

    """
    Time is the same as V times the number of inserts which is O(logV) so this function is O(VlogV) to create
    Space is O(V) since it adds n items to a list.
    """
    def __init__(self, nodes):
        self.min_index = 0
        self.heap = []
        self.values = []
        self.index_map = list(nodes.keys())
        # self.heapify(nodes)
        for node_id, value in nodes.items():
            self.insert(node_id, value)

    """
    O(1) time and space - helper function
    """
    def size(self):
        return len(self.heap)

    """
    Time is O(1) for inserts into lists, but O(logV) for bubble up so this is O(logV)
    Space is O(1) since all of these are O(1) space operations
    """
    def insert(self, node, value):
        self.heap.append(node)
        self.values.append(value)
        self.index_map[node] = self.size() - 1
        self.bubble_up(self.size() - 1)

    """
    Bubble up is at worst O(logV) time to go all the way up.  The rest are constant in this function
    Space is constant since it just keeps a constant number of values
    """

    def bubble_up(self, last_node_num):
        still_bubbling = True
        while still_bubbling and last_node_num != 0:
            if self.values[last_node_num] < self.get_parent_values(last_node_num):
                parent_index = self.get_parent_index(last_node_num)
                self.swap(last_node_num, parent_index)
                still_bubbling = True
                last_node_num = parent_index
            else:
                still_bubbling = False
        return

    """
    Time is O(logV) since bubble_up is O(logV) and the rest are constant
    Space is O(1), just look ups
    """
    def decrease_key(self, node_id, value):
        # find and decrease node.  O(logn) to find node
        index_of_node = self.index_map[node_id]
        if index_of_node is None:
            # node has already been visited -> no need to update it in the heap
            return
        self.values[index_of_node] = value
        # percolate up
        self.bubble_up(index_of_node)

    """
    All of the operations are O(1) to access and delete but sift down is O(logV) so this function is too.
    Space is O(1) since all memory usages are constant
    """
    def delete_min(self):
        min = self.heap[self.min_index]
        # base case - one left
        if len(self.heap) == 1:
            del self.heap[0]
            return min

        # update map
        self.index_map[self.heap[-1]] = 0
        self.index_map[self.heap[0]] = None

        # put last element first
        last_element = self.heap[-1]
        last_element_value = self.values[-1]
        self.heap[self.min_index] = last_element
        self.values[self.min_index] = last_element_value
        del self.heap[-1]
        del self.values[-1]

        # fix heap
        self.sift_down(self.min_index)
        return min

    """
    It is O(1) time since all it does it access and compare which are O(1) operations
    O(1) space since it only stores a few values
    """
    def get_min_child_value_and_index(self, node_num):
        left_child = float("inf")
        right_child = float("inf")
        if ((2 * node_num) + 1) < len(self.values):
            left_child = self.values[(2 * node_num) + 1]
        if ((2 * node_num) + 2) < len(self.values):
            right_child = self.values[(2 * node_num) + 2]
        if left_child <= right_child:
            return left_child, (2 * node_num) + 1
        else:
            return right_child, (2 * node_num) + 2

    """
    helper function O(1) time and space - one calculation and one look up
    """
    def get_parent_values(self, node_number):
        return self.values[(node_number - 1) // 2]

    """
    helper function O(1) time and space - one calculation
    """
    def get_parent_index(self, node_number):
        return (node_number - 1) // 2

    """
    This function is O(logv) worst case to sift a node all the way down the logv levels
    O(1) space since it only stores a constant amount of values.
    """
    def sift_down(self, node_num):
        while node_num < (len(self.heap) // 2):
            min_child, child_index = self.get_min_child_value_and_index(node_num)
            if min_child < self.values[node_num]:
                self.swap(child_index, node_num)
                # follow the node down and keep sifting
                node_num = child_index
            else:
                return

    """
    Swap is O(1) time since all it does it access array
    it is O(1) space since it just swaps them
    """
    def swap(self, child_index, parent_index):
        # all of these operations are O(1) just changing numbers

        # swap places in map
        parent_node = self.heap[parent_index]
        chid_node = self.heap[child_index]
        self.index_map[parent_node] = child_index
        self.index_map[chid_node] = parent_index

        # swap ids
        temp = self.heap[child_index]
        self.heap[child_index] = self.heap[parent_index]
        self.heap[parent_index] = temp

        # swap values
        temp = self.values[child_index]
        self.values[child_index] = self.values[parent_index]
        self.values[parent_index] = temp




