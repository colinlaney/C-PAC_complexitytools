
"""
This module implements th infomap community detection method
"""

__all__     = ["get_best_module", "module_hierachy"]
__authors__ = """Florian Gesser (gesser.florian@googlemail.com)"""

import math
import random
import sys

import numpy as np
import networkx as nx

from itertools import groupby


EXIT            = 'EXIT'
EPSILON_REDUCED = 1.0e-10
PASS_MAX        = sys.maxint #2^63 - 1 on 64bit machines

"""
Class that abstracts and capsules the Infomap auxiliary datastructures 
and provides the method interface to compute the Infomap algorithm on an undriected input graph.
"""
class Infomap(object):
    final_codelength = 0.0
    """ 
    Class initializer to setup proper standard values for the used datastructures
    """
    def __init__(self, graph):
        super(Infomap, self).__init__()
        self.graph         = graph

        #pruge self loops
        loop_edges = self.graph.selfloop_edges()
        self.graph.remove_edges_from(loop_edges)

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


    def init(self):
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

        """ Exit is a proxy for the degree of the node, includding weights on self loops
            Use degrees as proxy for node frequency which is obtained from markov stationary distribution (random surfer calculation via page rank)

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


        self.exit = self.plogp(self.exitDegree)
        self.code_length = self.exit - 2.0 * self.exit_log_exit + self.degree_log_degree - self.nodeDegree_log_nodeDegree


    """
    This method gets called everytime after the second pass has run.
    It ensures that the bookkeeping data structures are properly set up, 
    as the Infomap algorithm computes on them

    Parameters
    ----------
    org_mods : dict
       a dictionary where the keys are the nodes and the values the associated modules

    Side effect
    -----------
    Computes the proper values of the infomap parameter varoables by taking the weight of the graph into account

    See Also
    --------
    __init__ and init do a similar setup. this method however takes into account the weight of the graph in its computation of the infomap parameters which is important after a new graph has been generated after each run of second phase.
    """
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


    """
    Computes the entropy using the inverse node degree
    in order to save division operations.
    This is a core method that gets called a lot and is therfore a critical
    path worth optimizing

    Parameters
    ----------
    degree : the empirical degree of the node which is being used in the calculation

    Returns
    -------
    The entropy H

    """
    def plogp(self, degree):
        """Entropy calculation"""
        p = self.inverseDegree * degree
        return p * math.log(p, 2) if degree > 0 else 0.0


    def get_random_permutation_of_nodes(self):
        nodes = self.graph.nodes()
        shuffeled_nodes = nodes[:]
        random.shuffle(shuffeled_nodes)
        return shuffeled_nodes


    """
    Computes the link strengh of the adjacent nodes and their
    corrosponding modules

    Parameters
    ----------
    node : The node whos neighbourhood should be investgated in therm of link
        in terms of link strenghts of the associated modules

    Returns
    -------
    weights : A dictionary with the computed weights for each neighbor module


    """
    def neighbourhood_link_strength(self, node):
        weights = {}
        for neighbor, datas in self.graph[node].items() :
            if neighbor != node :
                weight = datas.get("weight", 1)
                neighborcom = self.modules[neighbor]
                weights[neighborcom] = weights.get(neighborcom, 0) + weight

        return weights


    """
    Does the renumbering of the modules dictionary in order to prepare it
    for the next iteration in the algorithm

    Parameters
    ----------
    current_modules : the dictionary of the current node-module associations

    Returns
    -------
    ret : the renumbered module dictionary

    See Aslo
    -------
    test_renumbering_modules in the test_infomap for an illustrative example

    """
    def renumber_modules(self, current_modules):
        ret = current_modules.copy()
        module_list       = current_modules.values()
        uniqe_module_list = sorted(set(module_list), key=lambda x: module_list.index(x))
        mapping           = dict(zip(uniqe_module_list, range(len(uniqe_module_list))))

        for key in current_modules.keys():
            ret[key] = mapping[current_modules[key]]

        return ret

    """
     Computes the final node->module association at the end of the algorithm 
     across the generated levels/hierarchy

     Parameters
     ----------
     module_hierarchy : a list of dictionaries
        Each element of the list represenets one level of the computed modules (node->module association for that granularity)

    level : an integer
        the number of levels in the hierarchy

    Returns
    -------
    Final mapping of nodes to modules

    """   
    def generate_module_mapping(self, module_hierarchy, level):
        module_mapping = module_hierarchy[0].copy()
        for index in range(1, level + 1):
            for node, module in module_mapping.items():
                module_mapping[node] = module_hierarchy[index][module]
        return module_mapping


    """
    The core of the greedy optimization step of the Infomap algorithm

    For each node this method investigates all neighbours of the node and
    calculates how much (if at all) the map eqation could be reduces if the
    node would be moved in that associated neighbours module.

    In the end the node gets moved to the neighbours module which gives the 
    minimal description length possible.

    See Also:
    ---------

    A detailed desciption of the procedure can of couese be found in the
    actual Infomap papers [1] as well as in my detailed blog post 
    "Theory, Mathematics, Implementation and Ground Truth" [2]

    [1] Rosvall, M., Axelsson, D., & Bergstrom, C. T. (2010). The map equation. The European Physical Journal Special Topics 178(1), 13-23
    [2] Theory, Mathematics, Implementation and Ground Truth, 
    http://fgesser.io/theory-mathematics-implementation-and-ground-truth.html#theory-mathematics-implementation-and-ground-truth on fgesser.io

    """
    def determine_best_new_module(self, iteration, stat):
        randomSequence = self.get_random_permutation_of_nodes()
        for index, curr_node in enumerate(self.graph):
            pick = curr_node
            # filter out self links
            Nlinks = len(list(filter(lambda x: x != pick, self.graph.neighbors(pick))))
            wNtoM  = self.neighbourhood_link_strength(pick)
            fromM  = self.modules[pick]
            #that would be wrong, it would sum up all the edges from the neighbour
            #wfromM = sum([self.graph.node[neighbour][EXIT] for neighbour in self.graph.neighbors(pick)])
            #instead what we want is to look up the weight to own module in the  dict containing the neighbours links
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


    """
    The "first psss" or "first phase" as it gets descripted in the
    papers.

    This phase encompases the iterated calling of "determine_best_new_module"
    till a certain epsilon is reached and the description length can't be 
    further reduced, wich in turn triggers the "second pass"


    See Also
    --------
    second_pass : the next step in the iteration


    """
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


    """
    The "second pass" or "second phase" as it gets described in the
    papers.

    Nodes belonging to the same community as detected in "first pass" are merged
    into a single node.
    The new graph is build up from these so called "super nodes"


    Parameters
    ----------

    current_partition : the old partition (node->module associations) from which
        the new graph will be build.

    Returns
    -------

    aggregated_graph : The newly aggregated graph consisting of the computed "super nodes".
        Edges between nodes who belong to the same module are turned into self loops of the
        corresponding suoer node 

    """
    def second_pass(self, current_partition):
        aggregated_graph = nx.Graph()

        # The new graph consists of as many "supernodes" as there are partitions
        aggregated_graph.add_nodes_from(set(current_partition.values()))
        # make edges between communites, bundle more edges between nodes in weight attribute
        edge_list=[(current_partition[node1], current_partition[node2], attr.get('weight', 1) ) for node1, node2, attr in self.graph.edges(data=True)]
        sorted_edge_list = sorted(edge_list)
        sum_z = lambda tuples: sum(t[2] for t in tuples)
        weighted_edge_list = [(k[0], k[1], sum_z(g)) for k, g in groupby(sorted_edge_list, lambda t: (t[0], t[1]))]
        aggregated_graph.add_weighted_edges_from(weighted_edge_list)

        return aggregated_graph

"""
Main iteration that calls all the other sub methods of that the Infomap algorithm consists.

The pattern is:

init
First Pass
Second Pass

Init
First Pass
Second Pass

Init
.
.
.

"""
def infomap_iteration(graph):
    iteration =0
    partition = Infomap(graph)
    partition.init()
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
    iteration += 1

    while True:
        partition.first_pass(iteration)
        new_codelength = partition.code_length
        Infomap.final_codelength = new_codelength
        if ( (codelength - new_codelength) < EPSILON_REDUCED) or (new_codelength < 1.0):
            break
        best_partition = partition.renumber_modules(partition.modules)
        parition_list.append(best_partition)
        codelength = new_codelength
        org_mods = partition.modules
        current_graph = partition.second_pass(best_partition)
        partition.graph = current_graph
        partition.set_up(org_mods)
        iteration += 1

    return parition_list[:], partition


def number_of_modules_detected(hierarchy):
    return len(set(hierarchy[-1].values()))

"""
Computes a hierarchy/dendrogram of the detected modules

Returns
-------

module_hierarchy : list of dictionaries

"""
def module_hierachy(graph):
    module_hierachy, handle = infomap_iteration(graph)
    return module_hierachy


"""
Computes the partition for which the Map Equation is minimal

Returns
-------

Dictionary with final node->module associations


"""
def get_best_module(graph):
    module_hierachy, handle = infomap_iteration(graph)
    return handle.generate_module_mapping(module_hierachy, len(module_hierachy)-1)

"""
Main entry method that gets called if Infomap gets called as a module
(e.g from the terminal)

Parameters
----------
graph : the graph in which the communities should be detected.

Returns
-------
- the result of the community detection, the best partition in modules for this graph
- the achived minimal codelength

"""
def main(graph=None):
    if graph == None:
        graph = nx.karate_club_graph()
        print "No input graph provided, running infomap with Zachary's Karate Club"

    partition   = Infomap(graph)
    print "Detecting modules by minimizing codelength"
    dendro      = module_hierachy(graph)
    best_module = get_best_module(graph)
    final = Infomap.final_codelength
    print "Done"
    print "Final codelength: " + str(final) + " in " + str(number_of_modules_detected(dendro)) + " Modules"

if __name__ == '__main__':
    main()