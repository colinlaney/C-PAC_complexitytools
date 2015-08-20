import unittest
import networkx as nx
import random

import infomap

def girvan_graphs(zout) :
    """
    Create a graph of 128 vertices, 4 communities, like in
    Community Structure in  social and biological networks.
    Girvan newman, 2002. PNAS June, vol 99 n 12

    Community is node modulo 4
    """

    pout = float(zout)/96.
    pin = (16.-pout*96.)/31.
    graph = nx.Graph()
    graph.add_nodes_from(range(128))
    for x in graph.nodes() :
        for y in graph.nodes() :
            if x < y :
                val = random.random()
                if x % 4 == y % 4 :
                    #nodes belong to the same community
                    if val < pin :
                        graph.add_edge(x, y)

                else :
                    if val < pout :
                        graph.add_edge(x, y)
    return graph


class InfomapTest(unittest.TestCase):

    def test_girvan(self):
        """
        Generate ground truth graph with 4 communities.
        Test that infomap detects them correctly
        """
        graph = girvan_graphs(4)
        modules = infomap.infomap(graph)

        for node, module in modules.items():
            self.assertEqual(module, modules[node%4])

    def test_communites_to_nodes(self):
        """
        Test that the new graph generated in the the second phase has the number of communities as nodes
        """
        graph = nx.erdos_renyi_graph(50, 0.1)
        partition = dict([])
        for node in graph.nodes() :
            partition[node] = node % 5
        #set up mock object
        info_partition = infomap.Partition(graph)
        info_partition.modules = partition

        self.assertSetEqual(set(partition.values()), set(info_partition.second_pass((info_partition.modules)).nodes()))

    def test_isolated_nodes(self) :
        """
        Test that the generated graph is isomoprh for the pathological case when all nodes are isolated
        """
        graph = nx.erdos_renyi_graph(50, 0.1)
        partition = dict([])
        for node in graph.nodes() :
            partition[node] = node

        #set up mock object
        info_partition = infomap.Partition(graph)
        info_partition.modules = partition
        generated_graph = info_partition.second_pass((info_partition.modules))

        self.assertTrue(nx.is_isomorphic(graph, generated_graph))

if __name__ == '__main__':
    unittest.main()