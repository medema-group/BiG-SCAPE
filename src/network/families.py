"""Contains functions to execute affinity propagaion on arbitrary graphs"""

# from python
import warnings

# from dependencies
import numpy as np
import networkx as nx
from networkx import Graph
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning

# from other modules
from src.data import DB
from src.comparison.binning import RecordPair

# from this module
from .utility import edge_list_to_adj_list, sim_matrix_from_adj_list


def generate_families(
    connected_component: list[tuple[int, int, float, float, float, float]],
) -> list[tuple[int, int, float]]:
    """Execute affinity propagation on a connected component

    Args:
        connected_component (list[tuple[int, int, float, float, float, float]]):
            connected component in the form of a list of edges

    Returns:
        list[tuple[int, int, float]]: list of (region_id, family, cutoff) tuples
    """
    adj_list = edge_list_to_adj_list(connected_component)

    # this list is going to be in the same order as the distance matrix rows/col
    # and the list of labels after AP
    node_ids = list(adj_list.keys())

    distance_matrix = sim_matrix_from_adj_list(adj_list)

    labels, centers = aff_sim_matrix(distance_matrix)

    # assemble list of (region_id, family, cutoff) tuples for easy insertion
    # into db
    regions_families = []

    for idx, label in enumerate(labels):
        label = int(label)
        if label == -1:
            continue

        region_id = node_ids[idx]

        family = int(label)

        regions_families.append((region_id, family, 0.3))

    return regions_families


def aff_sim_matrix(matrix):
    """Execute affinity propagation on a __similarity__ matrix

    Note: a similarity matrix. Not a distance matrix.

    Args:
        matrix (numpy.array[numpy.array]): similarity matrix in numpy array of array
        format.

    Returns:
        tuple[list[int], list[int]]: list of labels and list of cluster center ids
    """
    # thanks numpy but we sort of know what we're doing
    warnings.filterwarnings(action="ignore", category=ConvergenceWarning)

    af = AffinityPropagation(
        damping=0.9,
        max_iter=1000,
        convergence_iter=200,
        affinity="precomputed",
    ).fit(matrix)

    return af.labels_, af.cluster_centers_indices_


def save_to_db(regions_families):
    """Save families to database

    Args:
        regions_families (list[tuple[int, int, float]]): list of (region_id, family,
        cutoff) tuples
    """
    bgc_record_family_table = DB.metadata.tables["bgc_record_family"]
    insert_statement = (
        bgc_record_family_table.insert()
        .values(regions_families)
        .prefix_with("OR REPLACE")
    )

    DB.execute(insert_statement)
