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

    def getShortestPath( self, destIndex ):
        self.dest = destIndex

        path_edges = []
        total_length = 0
        nodes = self.network.nodes
        prev = self.node_info["prev"]
        cost = self.node_info["cost"]
        current_node_index = destIndex
        while current_node_index != self.source:
            prev_node = prev[current_node_index]
            if prev_node is None:
                print("There is no previous node.  Quitting")
                return {'cost': float("inf"), 'path': []}
            edge_length = self.get_node_length(prev_node, current_node_index)
            if edge_length is None:
                print("There is no edge.  Quitting")
                return {'cost': float("inf"), 'path': []}
            path_edges.append((nodes[current_node_index].loc, nodes[prev_node].loc, '{:.0f}'.format(edge_length)))
            total_length += edge_length
            current_node_index = prev_node
        assert (int(total_length) == int(cost[destIndex]))
        return {'cost': cost[destIndex], 'path': path_edges}

    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        self.node_info = self.dijkstras(self.network.nodes, self.source, use_heap)
        t2 = time.time()
        return (t2-t1)

    def get_node_length(self, current_node_index, prev_node):
        for edge in self.network.nodes[current_node_index].neighbors:
            if edge.dest.node_id == prev_node:
                return edge.length

    def dijkstras(self, graph, start_node, use_heap):
        try:
            # initialize
            dist = {}
            prev = {}
            length_dict = {}
            for index, node in enumerate(graph):
                dist[node.node_id] = float('inf')
                prev[node.node_id] = None

            # run dijkstras
            dist[start_node] = 0
            if use_heap:
                priority = MinHeap(dist)
                if priority.index_map[start_node] != 0:
                    print("Heaping is not working")
            else:
                priority = Queue(start_node, dist.keys())

            iteration = 0
            # Checking size is O(n) worst case for list implementation and O(1) for heap
            while priority.size() != 0:
                # delete min is
                u = priority.delete_min()
                for edge in graph[u].neighbors:
                    if edge.length + dist[u] < dist[edge.dest.node_id]:
                        dist[edge.dest.node_id] = dist[u] + edge.length
                        prev[edge.dest.node_id] = u
                        priority.decrease_key(edge.dest.node_id, dist[edge.dest.node_id])
                iteration += 1
            return {"cost": dist, "prev": prev}

        except Exception as e:
            print(e)
            traceback.print_tb(e.__traceback__)

# Heap needs to have access to array - maybe have two arrays to see where it is in the queue

class Queue:

    """
    creating the full list is n times of inserting, which is O(n)
    Space is O(n) since it stores all the nodes
    """
    def __init__(self, start_node, list_of_all_nodes):
        # self.queue = [float("inf")] * len(list_of_all_nodes)
        self.queue = []
        # insert for loop -> n times O(1) or O(n)
        for node_id in list_of_all_nodes:
            self.insert(node_id, float("inf"))

        # for this lab set the start_node to zero distance
        self.queue[start_node] = 0

    """
    O(n) time and space worst case to go through each x in the queue and adding it to another list
    Note: this is never used in this data structure -> only outside
    """
    def size(self):
        return len([x for x in self.queue if x is not None])

    def delete_min(self):
        # this is O(n) since it takes n to find the min and n to find the index
        min_index = self.queue.index(min(x for x in self.queue if x is not None))
        min_node_num = self.queue[min_index]
        self.queue[min_index] = None
        return int(min_index)

    def decrease_key(self, key, value):
        # changing the value is O(1) for lookups
        self.queue[key] = value
        return

    def insert(self, key, value):
        self.queue.append(key)
        self.queue[key] = value
        return

class MinHeap:

    """
    Time is the same as n times the number of inserts which is O(logn) so this function is O(nlogn) to create
    Space is O(n) since it adds n items to a list.
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
    Time is O(1) for inserts into lists, but O(logn) for bubble up so this is O(logn)
    Space is O(1) since all of these are O(1) space operations
    """
    def insert(self, node, value):
        self.heap.append(node)
        self.values.append(value)
        self.index_map[node] = self.size() - 1
        self.bubble_up(self.size() - 1)

    """
    Bubble up is at worst O(logn) time to go all the way up.  The rest are constant in this function
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
    debug function O(n) but not used in main code now
    """

    def print(self):
        repeat_num = 1
        repeat_again = 1
        print("Start Heap")
        for index in range(len(self.heap)):
            if index == repeat_num:
                print("\n", end="")
                print("Node {}: value {} ,".format(self.heap[index], self.values[index]), end="")
                repeat_again *= 2
                repeat_num += repeat_again
            else:
                print("Node {}: value {} ,".format(self.heap[index], self.values[index]), end="")
        print("\n End Heap")

    """
    Time is O(logn) since bubble_up is O(logn) and the rest are constant
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
    All of the operations are O(1) to access and delete but sift down is O(logn) so this function is too.
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
    This function is O(1) worst case TODO: ask
    Space complexity is O(logn) to sift all the way down
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
    This function calls sift down log(n) times which means to create the heap is O(log^2(n))
    space is O(1) since it is constant in this and in sift down
    """

    def heapify(self, nodes):
        # only runs log(n) times for each layer
        i = len(nodes) // 2
        self.heap = list(nodes.keys())
        self.values = list(nodes.values())
        while i >= 0:
            self.sift_down(i)
            i -= 1

    # call insert function TODO: fix

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




