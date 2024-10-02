![License](https://img.shields.io/github/license/medema-group/BiG-SCAPE)
![Github downloads](https://img.shields.io/github/downloads/medema-group/BiG-SCAPE/latest/total?label=Github%20downloads%20%28latest%29)
![Conda downloads](https://img.shields.io/conda/dn/bioconda/bigscape?label=Conda%20downloads)
![Test workflow](https://github.com/medema-group/BiG-SCAPE/actions/workflows/test.yml/badge.svg)
![Coverage](https://medema-group.github.io/BiG-SCAPE/badges/coverage.svg)
![Docs workflow](https://github.com/medema-group/BiG-SCAPE/actions/workflows/docs.yml/badge.svg)
![Pylint](https://medema-group.github.io/BiG-SCAPE/badges/pylint.svg)


## Running BiG-SCAPE

There are a few ways to run BiG-SCAPE, depending on your needs.

### Prerequisites

Software:

- python 3.11 or up
- conda/mamba

All steps start with cloning the repository:

`git clone https://github.com/medema-group/BiG-SCAPE`.

All other steps assume you are in the BiG-SCAPE directory, using cd `BiG-SCAPE`.

The most straightforward ways of installing are using Docker or Hatch.
However, these do not allow you to run the application from anywhere using a command like `big-scape`.
If you want to do that, follow the **Using Conda/Mamba** steps.

### Using Docker

TODO

### Using Hatch

Install Hatch and setup an environment:

1. `pip install hatch`
2. `hatch -v env create bigscape`

Run the application through hatch in the BiG-SCAPE folder:

3. `hatch run bigscape:bigscape --help`

This command only works inside the BiG-SCAPE folder.


### Using Conda/Mamba

These steps use `mamba`.

1. `mamba env create -f environment.yml`
2. `mamba activate bigscape`
3. `pip install .`

You can now run bigscape anywhere:

`bigscape --help`


## Dev instructions

1. pip install -r dev-requirements.txt
2. pre-commit install
3. mkdocs build
4. have fun!
