"""Module containing a class for a BiG-SCAPE input object, which has all the
input parameters/arguments"""

# from python
import logging
from typing import Optional
from pathlib import Path

# from other modules
from src.errors import InvalidInputArgError


class InputParameters:
    """
    Class to store all run input parameters

    Attributes:
        gbk_path: Path
        input_mode: str
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
        self.input_mode: Optional[str] = None
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

    def parse(
        self,
        gbk_path: Path,
        input_mode: str,
        dataset_path: Path,
        metadata_path: Path,
        pfam_path: Path,
        download_pfam: bool,
        mibig_version: str,
        reference_path: Path,
        include_gbk: str,
        exclude_gbk: str,
        query_bgc_path: Path,
        min_bgc_size: int,
    ):
        """Load input arguments from commandline ArgParser object

        Args:
            gbk_path (Path): Path to input gbk folder
            readin_mode (str): where to look for input gbk files
            dataset_path (Path): Path to dataset description tsv
            metadata_path (Path): Path to metadata tsv
            pfam_path (Path): Path to pressed Pfam files
            download_pfam (Bool): Toggle to download Pfam
            mibig_version (str): MiBIG version number
            reference_path (Path): Path to user provided reference gbks
            include_gbk (str): include only gbks with this string in name
            exclude_gbk (str): exclude gbks with this string in name
            query_bgc_path (Path): Path to query BGC gbk
            min_bgc_size (int): mininum length of a BGC to analyse
        """

        if not gbk_path.exists():
            logging.error("Path to input GBK dir is not valid")
            raise InvalidInputArgError()
        self.gbk_path = gbk_path
        # self.mode = mode

        if dataset_path and not dataset_path.exists():
            logging.error("Path to dataset tsv file is not valid")
            raise InvalidInputArgError()
        self.dataset_path = dataset_path

        if metadata_path and not metadata_path.exists():
            logging.error("Path to metadata tsv file is not valid")
            raise InvalidInputArgError()
        self.metadata_path = metadata_path

        if not pfam_path.exists():
            logging.error("Path to Pfam dir is not valid")
            raise InvalidInputArgError()
        self.pfam_path = pfam_path

        if reference_path and not reference_path.exists():
            logging.error("Path to reference GBK dir is not valid")
            raise InvalidInputArgError()
        self.reference_path = reference_path

        if query_bgc_path and not query_bgc_path.exists():
            logging.error("Path to query BGC GBK file is not valid")
            raise InvalidInputArgError()
        self.query_bgc_path = query_bgc_path

        self.input_mode = input_mode
        self.download_pfam = download_pfam
        self.mibig_version = mibig_version
        self.include_gbk = include_gbk
        self.exclude_gbk = exclude_gbk
        self.mib_bgc_size = min_bgc_size
