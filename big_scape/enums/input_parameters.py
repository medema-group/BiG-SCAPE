"""Contains enums for input parameters."""

# from python
from enum import Enum


class INPUT_MODE(Enum):
    FLAT = "flat"
    RECURSIVE = "recursive"


class RUN_MODE(Enum):
    CLUSTER = "cluster"
    QUERY = "query"
    BENCHMARK = "benchmark"
    DEREPLICATE = "dereplicate"
