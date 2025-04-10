"""Contains code relevant to the GBK component"""

# from python
from __future__ import annotations
from pathlib import Path

# from dependencies

# from other modules
import big_scape.enums as bs_enums

# from this module


class GBK:
    """
    Class to describe a GBK component

    Attributes:
        path (Path): path to the GBK file
        hash (str): content based hash of the entire GBK file
        components (Dict): dictionary of components in the GBK file
    """

    def __init__(self, path, hash, nt_length, as_version, source_type) -> None:
        self.name: bs_enums.COMPONENTS = bs_enums.COMPONENTS.GBK
        self.path: Path = path
        self.hash: str = hash
        self.nt_length: int = nt_length
        self.as_version: str = as_version
        self.source_type: bs_enums.SOURCE_TYPE = source_type
        self.components: dict[str, list] = {}

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
