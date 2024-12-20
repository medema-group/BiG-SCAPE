"""Utility functions for network analysis"""

# from dependencies
import numpy as np


def edge_list_to_adj_list(
    edge_list: list[tuple[int, int, float, float, float, float, int]]
) -> dict[int, dict[int, float]]:
    """Return an adjacency list from a list of edges

    values in the adjacency list are the similarity between the two regions, not
    the distance

    Args:
        edge_list (list[tuple[int, int, float, float, float, float, str]]): list of edges

    Returns:
        dict[int, dict[int, float]]: adjacency list
    """
    # set up the dictionary of dictionaries
    adj_list: dict[int, dict[int, float]] = {}
    for edge in edge_list:
        adj_list[edge[0]] = {}
        adj_list[edge[1]] = {}

    # go through all edges
    for record_a, record_b, distance, _, _, _, _ in edge_list:
        adj_list[record_a][record_b] = 1 - distance
        adj_list[record_b][record_a] = 1 - distance

    # add any missing adjacencies by setting them to similarity of 0
    for record_a in adj_list:
        for record_b in adj_list:
            if record_b not in adj_list[record_a]:
                adj_list[record_a][record_b] = 0.0
            if record_a not in adj_list[record_b]:
                adj_list[record_b][record_a] = 0.0

        # self similarity
        adj_list[record_a][record_a] = 1.0

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
    region_to_index = {}
    for index, region in enumerate(adj_list):
        region_to_index[region] = index

    # go through all regions
    for record_a in adj_list:
        a_matrix_idx = region_to_index[record_a]
        for record_b in adj_list[record_a]:
            b_matrix_idx = region_to_index[record_b]
            matrix[a_matrix_idx][b_matrix_idx] = 1 - adj_list[record_a][record_b]

    return matrix.tolist()


def edge_list_to_sim_matrix(
    edge_list: list[tuple[int, int, float, float, float, float, int]]
) -> tuple[np.ndarray, list[int]]:
    """Return a similarity matrix and id list from an edge list

    This function returns an adjacency matrix, and a list where each element is the id
    of a node at that index in the matrix

    Args:
        edge_list (list[tuple[int, int, float, float, float, float, str]]): list of edges

    Returns:
        tuple[np.ndarray, list[int]]: similarity matrix and list of region ids
    """

    # making sure this is correct by going through the edge list twice
    # first time to get all the ids
    id_to_idx: dict[int, int] = {}
    idx_to_id: list[int] = []

    for a_idx, b_idx, _, _, _, _, _ in edge_list:
        if a_idx not in id_to_idx:
            id_to_idx[a_idx] = len(id_to_idx)
            idx_to_id.append(a_idx)
        if b_idx not in id_to_idx:
            id_to_idx[b_idx] = len(id_to_idx)
            idx_to_id.append(b_idx)

    # second time to fill in the matrix
    matrix = np.zeros((len(id_to_idx), len(id_to_idx)))

    for a_idx, b_idx, distance, _, _, _, _ in edge_list:
        a_matrix_idx = id_to_idx[a_idx]
        b_matrix_idx = id_to_idx[b_idx]
        matrix[a_matrix_idx][b_matrix_idx] = 1 - distance
        matrix[b_matrix_idx][a_matrix_idx] = 1 - distance

    # self similarities
    for i in range(len(matrix)):
        matrix[i][i] = 1.0

    return matrix, idx_to_id
