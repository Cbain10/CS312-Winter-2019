#!/usr/bin/python3


from CS312Graph import *
import time





import traceback

class NetworkRoutingSolver:
    def __init__( self):
        pass

    def get_node_length(self, current_node_index, prev_node):
        for edge in self.network.nodes[current_node_index].neighbors:
            if edge.dest.node_id == prev_node:
                return edge.length

    def initializeNetwork( self, network ):
        assert( type(network) == CS312Graph )
        self.network = network

    def getShortestPath( self, destIndex ):
        self.dest = destIndex
        # TODO: RETURN THE SHORTEST PATH FOR destIndex
        #       INSTEAD OF THE DUMMY SET OF EDGES BELOW
        #       IT'S JUST AN EXAMPLE OF THE FORMAT YOU'LL 
        #       NEED TO USE
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
                return None
            edge_length = self.get_node_length(prev_node, current_node_index)
            if edge_length is None:
                print("There is no edge.  Quitting")
                return None
            path_edges.append((nodes[current_node_index].loc, nodes[prev_node].loc,  '{:.0f}'.format(edge_length)))
            total_length += edge_length
            current_node_index = prev_node
        assert(int(total_length) == int(cost[destIndex]))
        return {'cost': cost[destIndex], 'path': path_edges}

    def computeShortestPaths( self, srcIndex, use_heap=False ):
        self.source = srcIndex
        t1 = time.time()
        self.node_info = self.dijkstras(self.network.nodes, self.source, use_heap)
        # TODO: RUN DIJKSTRA'S TO DETERMINE SHORTEST PATHS.
        #       ALSO, STORE THE RESULTS FOR THE SUBSEQUENT
        #       CALL TO getShortestPath(dest_index)
        t2 = time.time()
        return (t2-t1)


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
            while priority.size() != 0:  # or Q not empty
                # u = queue.pop()
                u = priority.delete_min()
                # TODO: insert?
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

    def __init__(self, start_node, list_of_all_nodes):
        self.queue = [float("inf")] * len(list_of_all_nodes)
        self.queue[start_node] = 0
        # self.index_map = list(list_of_all_nodes)

    def size(self):
        return len([x for x in self.queue if x is not None])

    def delete_min(self):
        # this is O(n) since it takes n to find the min and n to find the index
        min_index = self.queue.index(min(x for x in self.queue if x is not None))
        min_node_num = self.queue[min_index]
        self.queue[min_index] = None
        # self.index_map[min_index] = None
        # # decrease the values to compensate
        # for value in self.index_map[min_index + 1:]:
        #     value = value - 1
        return int(min_index)

    def decrease_key(self, key, value):
        # changing the value is O(1) for lookups
        # real_key = self.index_map[key]
        self.queue[key] = value
        return





class MinHeap:

    def __init__(self, nodes):
        self.min_index = 0
        self.heap = []
        self.values = []
        self.index_map = list(nodes.keys())
        self.heapify(nodes)

    def size(self):
        return len(self.heap)

    def insert(self, node):
        self.heap.append(node.id)
        self.bubble_up(self.size())

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


    def decrease_key(self, node_id, value):
        # find and decrease node.  O(logn) to find node
        index_of_node = self.index_map[node_id]
        if index_of_node is None:
            # node has already been visited -> no need to update it in the heap
            return
        self.values[index_of_node] = value
        # percolate up
        self.bubble_up(index_of_node)


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

    def get_min_child_value_and_index(self, node_num):
        left_child = float("inf")
        right_child = float("inf")
        if ((2*node_num)+1) < len(self.values):
            left_child = self.values[(2*node_num)+1]
        if ((2 * node_num) + 2) < len(self.values):
            right_child = self.values[(2*node_num)+2]
        if left_child <= right_child:
            return left_child, (2*node_num)+1
        else:
            return right_child, (2*node_num)+2

    def get_parent_values(self, node_number):
        return self.values[(node_number-1) // 2]

    def get_parent_index(self, node_number):
        return (node_number-1) // 2

    def sift_down(self, node_num):
        while node_num < (len(self.heap) // 2):
            min_child, child_index = self.get_min_child_value_and_index(node_num)
            if min_child < self.values[node_num]:
                self.swap(child_index, node_num)
                # follow the node down and keep sifting
                node_num = child_index
            else:
                return

    def heapify(self, nodes):
        # since we can leave the leaves alone, start on the top half
        i = len(nodes) // 2
        self.heap = list(nodes.keys())
        self.values = list(nodes.values())
        while i >= 0:
            self.sift_down(i)
            i -= 1

    def swap(self, child_index, parent_index):
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


