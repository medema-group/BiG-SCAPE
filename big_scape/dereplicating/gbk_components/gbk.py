"""Contains code relevant to the GBK component"""

# from python
from __future__ import annotations
from pathlib import Path
import logging

# from dependencies

# from other modules
import big_scape.enums as bs_enums
from big_scape.cli.config import BigscapeConfig

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
        self.enum: bs_enums.COMPONENTS = bs_enums.COMPONENTS.GBK
        self.path: Path = path
        self.hash: str = hash
        self.nt_length: int = nt_length
        self.as_version: str = as_version
        self.source_type: bs_enums.SOURCE_TYPE = source_type
        self.components: dict[bs_enums.COMPONENTS, list] = {}

    @staticmethod
    def length_filter(gbk_list: list[GBK]) -> list[GBK]:
        """Filters out GBK objects that are too long or too short

        Args:
            gbk_list (list[GBK]): list of GBK objects

        Returns:
            list[GBK]: filtered list of GBK objects
        """

        orig_size = len(gbk_list)
        filtered_gbks = [
            gbk
            for gbk in gbk_list
            if BigscapeConfig.MIN_BGC_LENGTH
            < gbk.nt_length
            < BigscapeConfig.MAX_BGC_LENGTH
        ]
        filtered_size = len(filtered_gbks)

        if filtered_size == 0:
            logging.error(
                "No GBKs remain after length constraint of min %s to max %s bp!",
                BigscapeConfig.MIN_BGC_LENGTH,
                BigscapeConfig.MAX_BGC_LENGTH,
            )
            raise RuntimeError()

        if orig_size != filtered_size:
            logging.info(
                "%s GBKs remain after length constraint of min %s to max %s bp",
                filtered_size,
                BigscapeConfig.MIN_BGC_LENGTH,
                BigscapeConfig.MAX_BGC_LENGTH,
            )
        return filtered_gbks

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
