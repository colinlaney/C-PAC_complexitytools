
"""
This module implements th infomap community detection method
"""
#__all__ = [""]
__author__ = """Florian Gesser (gesser.florian@googlemail.com)"""

NODE_FREQUENCY  = 'NODE_FREQUENCY'
EXIT            = 'EXIT'
EPSILON_REDUCED = 1.0e-10
PASS_MAX        = -1

import math
from itertools import groupby
import sys

import numpy as np
import networkx as nx

sys.path.append("/Users/florian/Data/Pending/GSOC/code/community_evaluation/mini_pipeline_community/")
import buildTestGraph as btg

class Partition(object):
    """Represents a partition of the graph"""
    def __init__(self, graph):
        super(Partition, self).__init__()
        self.graph         = graph
        #self.modules = dict(zip(self.graph, range(self.graph.nodes()[0], graph.number_of_nodes())))
        self.modules = dict([])

        count = 0
        for node in self.graph.nodes():
            self.modules[node] = count
            count += 1

        self.Nnode         = self.graph.number_of_nodes()
        self.Nmod          = self.Nnode
        self.degree        = sum(self.graph.degree().values())
        self.inverseDegree = 1.0/self.degree


        self.nodeDegree_log_nodeDegree = 0.0
        self.exit_log_exit     = 0.0
        self.degree_log_degree = 0.0
        self.exitDegree        = 0.0
        self.exit 			   = 0.0
        self.code_length       = 0.0

        self.mod_exit    = dict([])
        self.mod_degree  = dict([])



    def init(self, at_beginning=False):

        #TODO  what has to change when this gets called in later iterations?

        self.modules = dict([])
        self.Nnode         = self.graph.number_of_nodes()
        self.Nmod          = self.Nnode
        self.degree        = sum(self.graph.degree().values())
        self.inverseDegree = 1.0/self.degree

        count = 0
        for node in self.graph.nodes():
            self.modules[node] = count
            count += 1

        """node frequency as computed by page rank currenlty not used"""
        page_ranks = nx.pagerank(self.graph)
        nx.set_node_attributes(self.graph, 'NODE_FREQUENCY', page_ranks)
        """ In the beginning exit := degree of the node
            Use degrees as proxy for node frequency which is obtained from markov stationary distribution (random surfer calculation via page rank)

            Later: Exit := totoal weight of links to other modules
        """
        degrees={i:self.graph.degree(i) for i in self.graph}
        nx.set_node_attributes(self.graph, 'EXIT', degrees)

        self.nodeDegree_log_nodeDegree = sum([self.plogp(self.graph.degree(node)) for node in self.graph])

        for index, node in enumerate(self.graph):
            node_i_exit   = self.graph.node[node][EXIT]
            node_i_degree = self.graph.degree(node)

            self.exit_log_exit     += self.plogp(node_i_exit)
            self.degree_log_degree += self.plogp(node_i_exit + node_i_degree)
            self.exitDegree        += node_i_exit

            self.mod_exit[index]    = node_i_exit
            self.mod_degree[index]  = node_i_degree

        if at_beginning == True:
            self.exit = self.plogp(self.exitDegree)
            self.code_length = self.exit - 2.0 * self.exit_log_exit + self.degree_log_degree - self.nodeDegree_log_nodeDegree

    def plogp(self, degree):
        """Entropy calculation"""
        p = self.inverseDegree * degree
        return p * math.log(p, 2) if degree > 0 else 0.0

    def get_random_permutation_of_nodes(self):
        nodes = self.graph.nodes()
        return np.random.permutation(nodes)


    def neighbourhood_link_strength(self, node):
        weights = {}
        for neighbor, datas in self.graph[node].items() :
            if neighbor != node :
                weight = datas.get("weight", 1)
                neighborcom = self.modules[neighbor]
                weights[neighborcom] = weights.get(neighborcom, 0) + weight

        return weights

    def renumber_modules(self, current_modules):
        ret = current_modules.copy()
        vals = set(current_modules.values())
        mapping = dict(zip(vals,range(len(vals))))

        for key in current_modules.keys():
            ret[key] = mapping[current_modules[key]]

        return ret

	for key in communities.keys():
		ret[key] = mapping[communities[key]]

	return ret

    def determine_best_new_module(self, iteration):
        randomSequence = self.get_random_permutation_of_nodes()

        modif = True
        nb_pass_done = 0

        curr_mod = self.code_length

        while modif and nb_pass_done != PASS_MAX:
            curr_mod = self.code_length
            modif = False
            nb_pass_done += 1

            for index, curr_node in enumerate(self.graph):
                #pick   = randomSequence[index]
                pick = curr_node
                Nlinks = len(self.graph.neighbors(pick))
                wNtoM  = self.neighbourhood_link_strength(pick)
                fromM  = self.modules[pick]
                #that is wrong, it would sum up all the edges from the neighbour
                #wfromM = sum([self.graph.node[neighbour][EXIT] for neighbour in self.graph.neighbors(pick)])
                #instead what we want is to look up the weight to own module in the community_links dict
                wfromM =  wNtoM.get(fromM, 0.0)

                bestM       = fromM
                best_weight = 0.0
                best_delta  = 0.0

                NmodLinks = len((wNtoM.keys()))

                for key, value in wNtoM.items():
                    toM  = key
                    wtoM = value

                    deltaL = 0

                    correction = 0

                    if toM != fromM:
                        node_i_exit   = self.graph.node[pick][EXIT]
                        node_i_degree = self.graph.degree(pick)

                        delta_exit = self.plogp(self.exitDegree - 2*wtoM + 2*wfromM) - self.exit

                        delta_exit_log_exit = - self.plogp(self.mod_exit[fromM + correction])                               \
                                              - self.plogp(self.mod_exit[toM + correction])                                 \
                                              + self.plogp(self.mod_exit[fromM + correction] - node_i_exit + 2*wfromM)      \
                                              + self.plogp(self.mod_exit[toM + correction] + node_i_exit - 2*wtoM)

                        delta_degree_log_degree = - self.plogp(self.mod_exit[fromM +correction ] + self.mod_degree[fromM +correction])                                          \
                                                  - self.plogp(self.mod_exit[toM + correction] + self.mod_degree[toM + correction])                                              \
                                                  + self.plogp(self.mod_exit[fromM +correction ] + self.mod_degree[fromM +correction] - node_i_exit - node_i_degree + 2*wfromM) \
                                                  + self.plogp(self.mod_exit[toM + correction] + self.mod_degree[toM + correction] + node_i_exit + node_i_degree - 2*wtoM)

                        deltaL = delta_exit - 2.0 * delta_exit_log_exit + delta_degree_log_degree

                    if deltaL < best_delta:
                        bestM = toM
                        best_weight = wtoM
                        best_delta = deltaL

                if bestM != fromM:
                    modif = True
                    node_i_exit   = self.graph.node[pick][EXIT]
                    node_i_degree = self.graph.degree(pick)


                    self.exitDegree        -= self.mod_exit[fromM + correction] + self.mod_exit[bestM + correction]
                    self.exit_log_exit     -= self.plogp(self.mod_exit[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction])
                    self.degree_log_degree -= self.plogp(self.mod_exit[fromM + correction] + self.mod_degree[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction] + self.mod_degree[bestM + correction])

                    self.mod_exit[fromM + correction]    -= node_i_exit - 2*wfromM
                    self.mod_degree[fromM + correction]  -= node_i_degree

                    self.mod_exit[bestM + correction]    += node_i_exit - 2*best_weight
                    self.mod_degree[bestM + correction]  += node_i_degree


                    self.exitDegree        += self.mod_exit[fromM + correction] + self.mod_exit[bestM + correction]
                    self.exit_log_exit     += self.plogp(self.mod_exit[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction])
                    self.degree_log_degree += self.plogp(self.mod_exit[fromM + correction] + self.mod_degree[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction] + self.mod_degree[bestM + correction])

                    self.exit = self.plogp(self.exitDegree)

                    self.code_length = self.exit - 2.0 * self.exit_log_exit + self.degree_log_degree - self.nodeDegree_log_nodeDegree

                    self.modules[pick] = bestM

            if (curr_mod - self.code_length ) < EPSILON_REDUCED:
                break

    def first_pass(self, iteration):
        #while passes_done != PASS_MAX
        self.determine_best_new_module(iteration)




    def second_pass(self):
        aggregated_graph = nx.Graph()

        # The new graph consists of as many "supernodes" as there are partitions
        aggregated_graph.add_nodes_from(set(self.modules.values()))
        # make edges between communites, bundle more edges between nodes in weight attribute
        edge_list=[(self.modules[node1], self.modules[node2], attr.get('weight', 1) ) for node1, node2, attr in self.graph.edges(data=True)]
        sorted_edge_list = sorted(edge_list)
        sum_z = lambda tuples: sum(t[2] for t in tuples)
        weighted_edge_list = [(k[0], k[1], sum_z(g)) for k, g in groupby(sorted_edge_list, lambda t: (t[0], t[1]))]
        aggregated_graph.add_weighted_edges_from(weighted_edge_list)

        return aggregated_graph



def infomap(graph):
    #import pdb; pdb.set_trace()

    # partition.move()

    iteration =0

    partition = Partition(graph)
    partition.init(True)

    parition_list = list()
    partition.first_pass(iteration)
    best_partition = partition.modules
    new_codelength = partition.code_length
    partition.modules = partition.renumber_modules(best_partition)
    parition_list.append(partition.modules)
    codelength = new_codelength
    current_graph = partition.second_pass()
    partition.graph = current_graph
    partition.init(True)

    iteration += 1

    while True:
        partition.first_pass(iteration)
        best_partition = partition.modules
        new_codelength = partition.code_length
        if codelength - new_codelength < EPSILON_REDUCED :
            break
        partition.modules = partition.renumber_modules(best_partition)
        parition_list.append(partition.modules)
        codelength = new_codelength
        current_graph = partition.second_pass()
        codelength = new_codelength
        partition.graph = current_graph
        partition.init(True)

        iteration += 1
    return parition_list[:]


def main():
    #test prep
    #graph = btg.build_graph()
    import networkx as nx
    graph = nx.karate_club_graph()
    # call to main algorithm method
    graph_partition = infomap(graph)
    print graph_partition
    print len(set(graph_partition[0].values()))

if __name__ == '__main__':
    main()