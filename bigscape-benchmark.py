#!/usr/bin/env python3

"""
BiG-SCAPE benchmark


Usage:   Please see `python bigscape-benchmark.py -h`

Example:
    python bigscape-benchmark.py -i ./bigscape_output -o ./results -g ./curated_GCFs.tsv
"""

from benchmark.main import run_benchmark

if __name__ == "__main__":
    run_benchmark()
