"""Module that contains helper functions specifically related to the bigscape version
"""

import tomllib
from importlib import metadata
from pathlib import Path


def get_bigscape_version() -> str:
    """Get the version of BiG-SCAPE.
    The way we retrieve the version is different depending on whether the package is
    installed or not.

    We need a dedicated library for this because the python community has not figured
    out that version numbers are pretty core to software development and there is no
    single place to put them. We want it to only be in the pyproject.toml file and not
    anywhere else, but this file is not available when installed as a package
    """
    # can we get to the pyproject.toml file?
    pyproject_toml = Path(__file__).parent.parent.parent / "pyproject.toml"

    if pyproject_toml.exists():
        with open(pyproject_toml, "rb") as fp:
            return tomllib.load(fp)["project"]["version"]

    # if not, we're probably running as a package. get the version of the currently
    # installed big-scape package

    return metadata.version("big-scape")
