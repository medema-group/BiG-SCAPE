"""Contains classes and methods to create and manipulate networks as used in BiG-SCAPE
"""
from .bs_network import BSNetwork
from .affinity_propagation import (
    sim_matrix_from_graph,
    sim_matrix_from_edge_list,
    aff_sim_matrix,
)

__all__ = [
    "BSNetwork",
    "sim_matrix_from_graph",
    "sim_matrix_from_edge_list",
    "aff_sim_matrix",
]
