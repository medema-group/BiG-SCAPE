"""Module containing code to load and store GBK files

Authors:
"""

from pathlib import Path


class GBK:
    """Class to describe a given GBK file"""

    path: Path
    metadata: dict[str, str]

    def __init__(self, path):
        self.path = path

    @staticmethod
    def parse_gbk(path: Path):
        """Parses a GBK file and returns a GBK object with all necessary information"""
        return GBK(path)
