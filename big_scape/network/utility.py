"""Utility functions for network analysis"""

# from dependencies
import numpy as np
import networkx as nx


def sim_matrix_from_graph(graph: nx.Graph, edge_property: str) -> np.ndarray:
    """Return a similarity matrix from a graph in the form of a numpy array

    Args:
        graph (Graph): graph
        edge_property (str): _description_

    Returns:
        ndarray: _description_
    """
    matrix = nx.to_numpy_array(graph, weight=edge_property, nonedge=1.0)
    # have to convert from distances to similarity
    matrix = 1 - matrix
    return matrix


def edge_list_to_adj_list(
    edge_list: list[tuple[int, int, float, float, float, float]]
) -> dict[int, dict[int, float]]:
    """Return an adjacency list from a list of edges

    values in the adjacency list are the similarity between the two regions, not
    the distance

    Args:
        edge_list (list[tuple[int, int, float, float, float, float]]): list of edges

    Returns:
        dict[int, dict[int, float]]: adjacency list
    """
    # set up the dictionary of dictionaries
    adj_list: dict[int, dict[int, float]] = dict()
    for edge in edge_list:
        adj_list[edge[0]] = dict()
        adj_list[edge[1]] = dict()

    # go through all edges
    for region_a, region_b, distance, _, _, _ in edge_list:
        adj_list[region_a][region_b] = 1 - distance
        adj_list[region_b][region_a] = 1 - distance

    # add any missing adjacencies by setting them to similarity of 0
    for region_a in adj_list:
        for region_b in adj_list:
            if region_b not in adj_list[region_a]:
                adj_list[region_a][region_b] = 0.0
            if region_a not in adj_list[region_b]:
                adj_list[region_b][region_a] = 0.0

        # self similarity
        adj_list[region_a][region_a] = 1.0

    return adj_list


def adj_list_to_sim_matrix(adj_list: dict[int, dict[int, float]]) -> np.ndarray:
    """Return a similarity matrix from an adjacency list

    Adjacency list is expected to be a dictionary of dictionaries, where the keys of
    the outer dictionary are region ids and the keys of the inner dictionaries are
    region ids that are adjacent to the outer key. The values of the inner
    dictionaries are the distances between the two regions

    This function will return a matrix where the columns and rows are ordered in the
    same order as the adjacency keys

    Args:
        adj_list (dict[int, dict[int, float]]): adjacency list

    Returns:
        np.ndarray: similarity matrix
    """
    # set up the matrix
    matrix = np.zeros((len(adj_list), len(adj_list)))

    # set up a dictionary to map region ids to matrix indices
    region_to_index = dict()
    for index, region in enumerate(adj_list):
        region_to_index[region] = index

    # go through all regions
    for region_a in adj_list:
        a_matrix_idx = region_to_index[region_a]
        for region_b in adj_list[region_a]:
            b_matrix_idx = region_to_index[region_b]
            matrix[a_matrix_idx][b_matrix_idx] = 1 - adj_list[region_a][region_b]

    return matrix.tolist()
