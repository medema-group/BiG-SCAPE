"""Utility functions for network analysis"""

# from dependencies
import numpy as np


def edge_list_to_sim_matrix(
    edge_list: list[tuple[int, int, float, float, float, float, int]],
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
