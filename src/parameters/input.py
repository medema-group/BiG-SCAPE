"""Module containing a class for a BiG-SCAPE input object, which has all the
input parameters/arguments"""

# from python
from typing import Optional
from pathlib import Path

# from dependencies
# from other modules
# from this module


class Input:
    """
    Class to store all run input parameters

    Attributes:
        gbk_path: Path
        mode: str
        dataset_path: Path
        metadata_path: Path
        pfam_path: Path
        download_pfam: Bool
        mibig_version: str
        reference_path: Path
        include_gbk: str
        exclude_gbk: str
        query_bgc_path: Path
        mib_bgc_size: int
    """

    def __init__(self) -> None:
        self.gbk_path: Optional[Path] = None
        self.mode: Optional[str] = None
        self.dataset_path: Optional[Path] = None
        self.metadata_path: Optional[Path] = None
        self.pfam_path: Optional[Path] = None
        self.download_pfam: Optional[bool] = None
        self.mibig_version: Optional[str] = None
        self.reference_path: Optional[Path] = None
        self.include_gbk: Optional[str] = None
        self.exclude_gbk: Optional[str] = None
        self.query_bgc_path: Optional[Path] = None
        self.mib_bgc_size: Optional[int] = None
