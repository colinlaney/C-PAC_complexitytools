import networkx as nx
import itertools

def build_graph():
	edges1, edges2 = build_edgelists()
	all_eges = edges1 + edges2
	graph = nx.from_edgelist(all_eges)
	return graph


def build_edgelists():
	first_cluster  = {1:[2, 3,5], 2:[1, 4, 5], 3:[1, 4, 5], 4:[2, 3, 5], 5:[1, 2, 3, 4]}
	second_cluster =  build_second_cluster(first_cluster)
	 #{6: [7, 8, 10], 7: [6, 9, 10], 8: [6, 9, 10], 9: [7, 8, 10], 10: [6, 7, 8, 9]}

	edgelist_first_cluster  = cartesianProduct(first_cluster)
	edgelist_second_cluster = cartesianProduct(second_cluster)

	outgoing_edge_1 = (2,8)
	outgoing_edge_2 = (8,2)

	edgelist_first_cluster.append(outgoing_edge_1)
	edgelist_second_cluster.append(outgoing_edge_2)

	return edgelist_first_cluster, edgelist_second_cluster


def build_second_cluster(node_dict):
	tuples = []
	for k, v in node_dict.items():
		tuples.append((k+5,[k+5 for k in v]))
	second = dict(tuples)
	return second


def cartesianProduct(node_dict):
	iter_edges = []
	for k, v in node_dict.items():
		iter_edges.append(itertools.product([k], v))
	return list(itertools.chain.from_iterable(iter_edges))


if __name__ == '__main__':
	build_graph()