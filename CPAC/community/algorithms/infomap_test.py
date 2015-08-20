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

    def test_isolated_nodes(self):
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

    def test_renumbering_modules(self):
        """
        Test that the renumbering of the module dictionary is correct
        """
        #input           = {1: 2, 2: 2, 3: 2, 4: 2, 5: 2, 6: 6, 7: 6, 8: 6, 9: 6, 10: 6}
        input           = {0: 17, 1: 17, 2: 12, 3: 12, 4: 10, 5: 16, 6: 16, 7: 12, 8: 32, 9: 12, 10: 10, 11: 17, 12: 12, 13: 12, 14: 32, 15: 32, 16: 16, 17: 17, 18: 32, 19: 17, 20: 32, 21: 17, 22: 32, 23: 25, 24: 27, 25: 25, 26: 29, 27: 27, 28: 31, 29: 29, 30: 32, 31: 31, 32: 32, 33: 32}
        #expected_output = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 1, 7: 1, 8: 1, 9: 1, 10: 1}
        expected_output = {0: 0, 1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 3, 7: 1, 8: 4, 9: 1, 10: 2, 11: 0, 12: 1, 13: 1, 14: 4, 15: 4, 16: 3, 17: 0, 18: 4, 19: 0, 20: 4, 21: 0, 22: 4, 23: 5, 24: 6, 25: 5, 26: 7, 27: 6, 28: 8, 29: 7, 30: 4, 31: 8, 32: 4, 33: 4}
        #set up mock object, we just need some graph here to create the instance
        graph = nx.karate_club_graph()
        info_partition = infomap.Partition(graph)
        result = info_partition.renumber_modules(input)

        self.assertDictEqual(result, expected_output)



if __name__ == '__main__':
    unittest.main()