
"""
This module implements th infomap community detection method
"""
#__all__ = [""]
__author__ = """Florian Gesser (gesser.florian@googlemail.com)"""



import math
from itertools import groupby
import sys

import numpy as np
import networkx as nx

sys.path.append("/Users/florian/Data/Pending/GSOC/code/community_evaluation/mini_pipeline_community/")
import buildTestGraph as btg


NODE_FREQUENCY  = 'NODE_FREQUENCY'
EXIT            = 'EXIT'
EPSILON_REDUCED = 1.0e-10
PASS_MAX        = sys.maxint #2^63 - 1 on 64bit machines


class Partition(object):
    """Represents a partition of the graph"""
    def __init__(self, graph):
        super(Partition, self).__init__()
        self.graph         = graph

        #pruge self loops
        loop_edges = self.graph.selfloop_edges()
        self.graph.remove_edges_from(loop_edges)

        #self.modules = dict(zip(self.graph, range(self.graph.nodes()[0], graph.number_of_nodes())))
        self.modules = dict([])

        count = 0
        for node in self.graph.nodes():
            self.modules[node] = count
            count += 1

        self.Nnode         = self.graph.number_of_nodes()
        self.Nmod          = self.Nnode
        self.degree        = sum(self.graph.degree(weight="weight").values())
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
        #pruge self loops
        loop_edges = self.graph.selfloop_edges()
        #self.graph.remove_edges_from(loop_edges)

        self.modules = dict([])
        self.Nnode         = self.graph.number_of_nodes()
        self.Nmod          = self.Nnode
        self.degree        = sum(self.graph.degree(weight="weight").values())
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
        degrees={i:self.graph.degree(i, weight="weight") for i in self.graph}
        nx.set_node_attributes(self.graph, 'EXIT', degrees)

        self.nodeDegree_log_nodeDegree = sum([self.plogp(self.graph.degree(node, weight="weight")) for node in self.graph])

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
        #nodes = self.graph.nodes()
        #return np.random.permutation(nodes)
        nodes = self.graph.nodes()
        shuffeled_nodes = nodes[:]
        import random
        random.shuffle(shuffeled_nodes)
        return shuffeled_nodes


    def neighbourhood_link_strength(self, node):
        weights = {}
        for neighbor, datas in self.graph[node].items() :
            if neighbor != node :
                weight = datas.get("weight", 1)
                neighborcom = self.modules[neighbor]
                weights[neighborcom] = weights.get(neighborcom, 0) + weight

        return weights

    def renumber_modules(self, current_modules):
        # ret = current_modules.copy()
        # vals = set(current_modules.values())
        # mapping = dict(zip(vals,range(len(vals))))
        #
        # for key in current_modules.keys():
        #     ret[key] = mapping[current_modules[key]]
        #
        # return ret
        #
        count = 0
        ret = current_modules.copy()
        new_values = dict([])

        for key in current_modules.keys() :
            value = current_modules[key]
            new_value = new_values.get(value, -1)
            if new_value == -1 :
                new_values[value] = count
                new_value = count
                count = count + 1
            ret[key] = new_value

        return ret



    def determine_best_new_module(self, iteration, stat):
        randomSequence = self.get_random_permutation_of_nodes()


        for index, curr_node in enumerate(self.graph):
        #for index, curr_node in enumerate(randomSequence):
            #pick   = randomSequence[index]
            pick = curr_node


            # if index == 0:
            #     pick = 31
            # elif index == 1:
            #     pick = 11
            # else:
            #     pick = randomSequence[index]

            # if index == 0:
            #     pick = 5
            # elif index == 1:
            #     pick = 1
            #pick = curr_node

            # TODO in second phase wrong link count due to self loops

            #Nlinks = len(self.graph.neighbors(pick))
            Nlinks = len(list(filter(lambda x: x != pick, self.graph.neighbors(pick))))

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


            # TODO randomize neighbor links + sort the nodes to exit mapping in second step?
            # zaehler = 0
            # for mod, weight in wNtoM.items():
            #     randPos = np.random.permutation(len(wNtoM)
            #      temp_mod = mod
            #      temp_weight = weight
            #      wNtoM[zaehler].

            import random
            # keys = wNtoM.keys()
            # random.shuffle(keys)
            # for key in keys:
            #     print wNtoM[key]

            # items = wNtoM.items()
            # random.shuffle(items)

            for key, value in wNtoM.items():
                toM  = key
                wtoM = value

                deltaL = 0

                correction = 0

                if toM != fromM:
                    node_i_exit   = self.graph.node[pick][EXIT]
                    #node_i_degree = self.graph.degree(pick)
                    node_i_degree = self.graph.degree(pick, weight="weight")

                    delta_exit = self.plogp(self.exitDegree - 2.*wtoM + 2.*wfromM) - self.exit

                    delta_exit_log_exit = - self.plogp(self.mod_exit[fromM + correction])                               \
                                          - self.plogp(self.mod_exit[toM + correction])                                 \
                                          + self.plogp(self.mod_exit[fromM + correction] - node_i_exit + 2.*wfromM)      \
                                          + self.plogp(self.mod_exit[toM + correction] + node_i_exit - 2.*wtoM)

                    delta_degree_log_degree = - self.plogp(self.mod_exit[fromM +correction ] + self.mod_degree[fromM +correction])                                          \
                                              - self.plogp(self.mod_exit[toM + correction] + self.mod_degree[toM + correction])                                              \
                                              + self.plogp(self.mod_exit[fromM +correction ] + self.mod_degree[fromM +correction] - node_i_exit - node_i_degree + 2.*wfromM) \
                                              + self.plogp(self.mod_exit[toM + correction] + self.mod_degree[toM + correction] + node_i_exit + node_i_degree - 2.*wtoM)

                    deltaL = delta_exit - 2.0 * delta_exit_log_exit + delta_degree_log_degree

                if deltaL < best_delta:
                    bestM = toM
                    best_weight = wtoM
                    best_delta = deltaL

            if bestM != fromM:
                node_i_exit   = self.graph.node[pick][EXIT]
                node_i_degree = self.graph.degree(pick)


                self.exitDegree        -= self.mod_exit[fromM + correction] + self.mod_exit[bestM + correction]
                self.exit_log_exit     -= self.plogp(self.mod_exit[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction])
                self.degree_log_degree -= self.plogp(self.mod_exit[fromM + correction] + self.mod_degree[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction] + self.mod_degree[bestM + correction])

                self.mod_exit[fromM + correction]    -= node_i_exit - 2.*wfromM
                self.mod_degree[fromM + correction]  -= node_i_degree

                self.mod_exit[bestM + correction]    += node_i_exit - 2.*best_weight
                self.mod_degree[bestM + correction]  += node_i_degree


                self.exitDegree        += self.mod_exit[fromM + correction] + self.mod_exit[bestM + correction]
                self.exit_log_exit     += self.plogp(self.mod_exit[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction])
                self.degree_log_degree += self.plogp(self.mod_exit[fromM + correction] + self.mod_degree[fromM + correction]) + self.plogp(self.mod_exit[bestM + correction] + self.mod_degree[bestM + correction])

                self.exit = self.plogp(self.exitDegree)

                self.code_length = self.exit - 2.0 * self.exit_log_exit + self.degree_log_degree - self.nodeDegree_log_nodeDegree

                self.modules[pick] = bestM
                stat['moved'] = True



    def first_pass(self, iteration):
        Nloops = 0
        stat = {'moved':True}
        while stat['moved'] == True:
            stat['moved'] = False
            old_codelength = self.code_length
            self.determine_best_new_module(iteration, stat)
            Nloops += 1
            if (old_codelength - self.code_length) < 1.0E-10:
                stat['moved'] = False




    def second_pass(self, current_partition):
        # aggregated_graph = nx.Graph()
        #
        # # The new graph consists of as many "supernodes" as there are partitions
        # aggregated_graph.add_nodes_from(set(self.modules.values()))
        # # make edges between communites, bundle more edges between nodes in weight attribute
        # edge_list=[(self.modules[node1], self.modules[node2], attr.get('weight', 1) ) for node1, node2, attr in self.graph.edges(data=True)]
        # sorted_edge_list = sorted(edge_list)
        # sum_z = lambda tuples: sum(t[2] for t in tuples)
        # weighted_edge_list = [(k[0], k[1], sum_z(g)) for k, g in groupby(sorted_edge_list, lambda t: (t[0], t[1]))]
        # aggregated_graph.add_weighted_edges_from(weighted_edge_list)
        #
        # return aggregated_graph

        ret = nx.Graph()
        ret.add_nodes_from(current_partition.values())

        for node1, node2, datas in self.graph.edges_iter(data = True) :
            weight = datas.get("weight", 1)
            com1 = current_partition[node1]
            com2 = current_partition[node2]
            w_prec = ret.get_edge_data(com1, com2, {"weight":0}).get("weight", 1)
            ret.add_edge(com1, com2, weight = w_prec + weight)
        return ret

    def set_up(self, org_mods):

        self.modules = dict([])

        count = 0
        for node in self.graph.nodes():
            self.modules[node] = count
            count += 1

        indicies = list(set(org_mods.values()))
        exit_degrees = {_:self.mod_exit[index] for _, index in enumerate(indicies)}

        nx.set_node_attributes(self.graph, 'EXIT', exit_degrees)

        self.mod_exit = dict([])
        self.mod_degree = dict([])

        self.exit_log_exit = 0.0
        self.degree_log_degree = 0.0
        self.exitDegree = 0.0

        for index, node in enumerate(self.graph):
            node_i_exit   = self.graph.node[node][EXIT]
            node_i_degree = self.graph.degree(node, weight="weight")

            self.exit_log_exit     += self.plogp(node_i_exit)
            self.degree_log_degree += self.plogp(node_i_exit + node_i_degree)
            self.exitDegree        += node_i_exit

            self.mod_exit[index]    = node_i_exit
            self.mod_degree[index]  = node_i_degree

            self.exit = self.plogp(self.exitDegree)
            self.code_length = self.exit - 2.0 * self.exit_log_exit + self.degree_log_degree - self.nodeDegree_log_nodeDegree

    def generate_module_mapping(self, module_hierarchy, level):
        module_mapping = module_hierarchy[0].copy()
        for index in range(1, level + 1):
            for node, module in module_mapping.items():
                module_mapping[node] = module_hierarchy[index][module]
        return module_mapping



def iteration_loop(graph):
    #import pdb; pdb.set_trace()

    # partition.move()

    iteration =0

    partition = Partition(graph)
    partition.init(True)

    #uncompressedCodeLength = partition.code_length

    parition_list = list()
    partition.first_pass(iteration)
    new_codelength = partition.code_length
    best_partition= partition.renumber_modules(partition.modules)
    parition_list.append(best_partition)
    codelength = new_codelength
    org_mods = partition.modules
    current_graph = partition.second_pass(best_partition)
    partition.graph = current_graph
    partition.set_up(org_mods)
    #partition.init(True)

    iteration += 1

    while True:
        partition.first_pass(iteration)
        new_codelength = partition.code_length
        if ( (codelength - new_codelength) < EPSILON_REDUCED) or (new_codelength < 1.0):
            ret_code = codelength
            break
        best_partition = partition.renumber_modules(partition.modules)
        parition_list.append(best_partition)
        codelength = new_codelength
        org_mods = partition.modules
        current_graph = partition.second_pass(best_partition)
        partition.graph = current_graph
        partition.set_up(org_mods)
        #partition.init(True)

        iteration += 1
    return parition_list[:], partition


def infomap(graph):
    module_hierachy, handle = iteration_loop(graph)

    return handle.generate_module_mapping(module_hierachy, len(module_hierachy)-1)

def main():
    #test prep
    #graph = btg.build_graph()
    import networkx as nx
    #graph = nx.karate_club_graph()
    #graph = nx.read_gpickle("/Users/florian/Desktop/testgraph/testgraph")

    #fh=open("/Users/florian/Desktop/com-amazon.ungraph.txt")
    #graph = nx.read_edgelist(fh, nodetype=int)

    graph = girvan(4)

    # call to main algorithm method
    mapping = infomap(graph)

    print mapping
    nx.set_node_attributes(graph, 'finalmodule', mapping)
    drawNetwork(graph)

    #print "Final Codelength: " + str(codelength)
    #print "Compressed by: " + str((100.0*(1.0-codelength/uncompressedCodeLength)))
    #print "Levels: " + str(len(graph_partition))
    #print "Modules found last level: " + str(len(set(graph_partition[len(graph_partition)-2].values())))
    #print "Modules found last level: " + str(len(set(graph_partition[len(graph_partition)-1].values())))

def girvan(zout):
    pout = float(zout)/96.
    pin = (16.-pout*96.)/31.
    graph = nx.Graph()
    graph.add_nodes_from(range(128))
    for x in graph.nodes() :
        for y in graph.nodes() :
            if x < y :
                val = np.random.random()
                if x % 4 == y % 4 :
                    #nodes belong to the same community
                    if val < pin :
                        graph.add_edge(x, y)

                else :
                    if val < pout :
                        graph.add_edge(x, y)
    return graph


def drawNetwork(G):
    import matplotlib.pyplot as plt

    import matplotlib.colors as colors
    # position map
    pos = nx.spring_layout(G)
    # community ids
    communities = [v for k,v in nx.get_node_attributes(G, 'finalmodule').items()]
    numCommunities = max(communities) + 1
    # color map from http://colorbrewer2.org/
    cmapLight = colors.ListedColormap(['#a6cee3', '#b2df8a', '#fb9a99', '#fdbf6f', '#cab2d6'], 'indexed', numCommunities)
    cmapDark = colors.ListedColormap(['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00', '#6a3d9a'], 'indexed', numCommunities)

    # edges
    nx.draw_networkx_edges(G, pos)

    # nodes
    nodeCollection = nx.draw_networkx_nodes(G,
        pos = pos,
        node_color = communities,
        cmap = cmapLight
    )
    # set node border color to the darker shade
    darkColors = [cmapDark(v) for v in communities]
    nodeCollection.set_edgecolor(darkColors)

    # Print node labels separately instead
    for n in G.nodes_iter():
        plt.annotate(n,
            xy = pos[n],
            textcoords = 'offset points',
            horizontalalignment = 'center',
            verticalalignment = 'center',
            xytext = [0, 2],
            color = cmapDark(communities[n])
        )

    plt.axis('off')
    # plt.savefig("karate.png")
    plt.show()

if __name__ == '__main__':
    main()