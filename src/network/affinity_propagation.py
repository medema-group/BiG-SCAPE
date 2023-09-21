"""Contains functions to execute affinity propagaion on arbitrary graphs"""

# from python
import warnings

# from dependencies
from numpy import ndarray
import networkx as nx
from networkx import Graph
from sklearn.cluster import AffinityPropagation
from sklearn.exceptions import ConvergenceWarning

# from other modules
from src.data import DB
from src.comparison.binning import RecordPair


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


def sim_matrix_from_graph(graph: Graph, edge_property: str) -> ndarray:
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


def sim_matrix_from_edge_list(
    edge_list: list[tuple[RecordPair, float, float, float, float]]
) -> ndarray:
    """Return a similarity matrix from a list of edges in the form of a numpy array

    Args:
        edge_list (list[tuple[BGCPair, float, float, float, float]]): list of edges

    Returns:
        ndarray: _description_
    """
    # get the number of regions
    num_regions = len(set([pair.region_a for pair in edge_list]))

    # initialize the matrix
    matrix = [[0.0 for _ in range(num_regions)] for _ in range(num_regions)]

    # populate the matrix
    for pair, distance, _, _, _ in edge_list:
        matrix[pair.region_a._db_id][pair.region_b._db_id] = 1 - distance
        matrix[pair.region_b._db_id][pair.region_a._db_id] = 1 - distance

    return matrix


def save_families_to_db(edge_list, cutoff, families) -> None:
    """Save families to the database

    Args:
        edge_list (list[tuple[BGCPair, float, float, float, float]]): list of edges
        labels (list[int]): list of labels
        centers (list[int]): list of centers
    """

    families_table = DB.metadata.tables["bgc_record_family"]

    # save the entry to the database
