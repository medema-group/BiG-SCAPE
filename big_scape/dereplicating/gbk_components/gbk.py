"""Contains code relevant to the GBK component"""

# from python
from __future__ import annotations
from pathlib import Path

# from dependencies

# from other modules

# from this module


class GBK:
    """
    Class to describe a GBK component

    Attributes:
        path (Path): path to the GBK file
        hash (str): content based hash of the entire GBK file
        components (Dict): dictionary of components in the GBK file
    """

    def __init__(self, path, hash) -> None:
        self.path: Path = path
        self.hash: str = hash
        self.components: dict = {}

    @classmethod
    def create(cls, path: Path, hash: str) -> GBK:
        """Creates a new GBK object and returns it with its path and hash

        Args:
            path (Path): path to the GBK file
            hash (str): content based hash of the entire GBK file

        Returns:
            GBK: new GBK object
        """

        return cls(path, hash)

    def __repr__(self) -> str:
        return f"{self.path}"

    def __hash__(self) -> int:
        return hash(self.hash)

    def __eq__(self, other) -> bool:
        if not isinstance(other, GBK):
            return False

        if self.hash is None or other.hash is None:
            return False

        return self.hash == other.hash
