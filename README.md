![License](https://img.shields.io/github/license/medema-group/BiG-SCAPE)
![Github downloads](https://img.shields.io/github/downloads/medema-group/BiG-SCAPE/latest/total?label=Github%20downloads%20%28latest%29)
![Conda downloads](https://img.shields.io/conda/dn/bioconda/bigscape?label=Conda%20downloads)
![Test workflow](https://github.com/medema-group/BiG-SCAPE/actions/workflows/test.yml/badge.svg)
![Coverage](https://medema-group.github.io/BiG-SCAPE/badges/coverage.svg)
![Docs workflow](https://github.com/medema-group/BiG-SCAPE/actions/workflows/docs.yml/badge.svg)
![Pylint](https://medema-group.github.io/BiG-SCAPE/badges/pylint.svg)


# BiG-SCAPE

**BiG-SCAPE** (Biosynthetic Gene Similarity Clustering and Prospecting Engine) is a software package, written in Python, that constructs sequence similarity networks of Biosynthetic Gene Clusters (BGCs) and groups them into Gene Cluster Families (GCFs). BiG-SCAPE does this by rapidly calculating a distance matrix between gene clusters based on a comparison of their protein domain content, order, copy number and sequence identity.

As input, BiG-SCAPE takes GenBank files from the output of [antiSMASH](https://antismash.secondarymetabolites.org) with BGC predictions, as well as reference BGCs from the [MIBiG repository](https://mibig.secondarymetabolites.org/). As output, BiG-SCAPE generates tab-delimited output files, as well as a rich HTML visualization.

In principle, BiG-SCAPE can also be used on any other gene clusters, such as pathogenicity islands, secretion system-encoding gene clusters, or even whole viral genomes.

If you find BiG-SCAPE useful, please [cite us]() [TODO V2].



## Running BiG-SCAPE

There are a few ways to run BiG-SCAPE, depending on your needs.

### Prerequisites

Software:

- python 3.11 or up
- conda/mamba

### Using Conda/Mamba

These steps use `mamba`.

Clone this repository:

1. `git clone https://github.com/medema-group/BiG-SCAPE`.
2. `cd BiG-SCAPE
3. `mamba env create -f environment.yml`
4. `mamba activate bigscape`
5. `pip install .`

You can now run bigscape anywhere:

`bigscape --help`

### Using Docker
Run BiG-SCAPE through docker using the docker run command:

```sh
docker run \
    --volume your_root_data_dir:/home/data \
    --detatch=false \
    --rm \
    --user=$(id -u):$(id -g) \
    TODO: DOCKER URL \
    # arguments from here are the same as using bigscape.py normally
    cluster \
    -i /home/data/your_input \
    -o /home/data/output_folder \
    -p /home/data/pfam_folder/Pfam-A.hmm
```

`your_root_data_dir` must be a parent folder of your input, your Pfam and where you want
to put your output.

For example:

```
/home/example/data
  ├ /input
  |    ├ experiment_1/
  |    |  └ sample_1/
  |    |      ├ a.region001.gbk
  |    |      └ a.region002.gbk
  |    └ unrelated_folder_do_not_use/
  ├ /output
  └ /pfam
      └ Pfam-A.hmm
```

Can use a command as such:

```
docker run \
    --volume /home/example/data:/home/data \
    --detatch=false \
    --rm \
    --user=$(id -u):$(id -g) \
    TODO: DOCKER URL \
    cluster \
    -i /home/data/input/experiment_1 \
    -o /home/data/output/experument_1 \
    -p /home/data/pfam/Pfam-A.hmm
```


## Dev instructions

We appreciate any contributions!
In order to setup a development environment, please follow the following steps:

### Using mamba

1. Install (micro)mamba
2. `git clone https://github.com/medema-group/BiG-SCAPE`.
3. `cd BiG-SCAPE`
4. `mamba env create -f environment.yml`
5. `mamba activate bigscape`
6. `pip install . --extra dev`
7. `pre-commit install`

To run BiG-SCAPE:
`python bigscape.py --help`

You can use the `bigscape --help` runnable directly, but this must be re-installed each
time you make a change to the codebase.

To run tests:
`python -m pytest`

When developing, make sure to use the mamba environment you created (called bigscape).
