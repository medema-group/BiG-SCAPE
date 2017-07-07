# Introduction

Bioinformatically, mining (meta)genomes for **Biosynthetic Gene Clusters** 
(BGCs) encoding specialized metabolites would entail identifying and annotating 
BGCs on the genome and taking additional steps to define a distance between BGCs
in order to map the BGC diversity in similarity networks. These similarity 
networks would graphically summarize the diversity of the BGCs, as well as 
contain multiple annotations to help identify novel compounds, make ecological 
correlations and so on.

## Defining a distance

BGCs are essentially a collection of genes that code for proteins that work 
together to produce a compound. These proteins are most likely the most 
important factor when it comes to the final structure of the compound. Thus, a 
good distance metric should use information on the similarity of the proteins 
between two BGCs.

In this project, three indices are combined to define a final distance metric 
between any given pair of BGCs:
* The Jaccard index (J): The ratio between the distinct shared and distinct 
unshared domain types between two BGCs.
* The DDS (Domain Duplicate Score): Measures the sequence similarity between the
domains of both BGCs, for each type of domain. When each BGC contains only one
copy of a certain domain, the sequence similarity can be obtained directly, 
otherwise, the Hungarian algorithm is used to select the most similar Pfam 
domain sequences (`munkres.py`). A special weight can be also be given to marked 
domains annotated as "anchor domains" (`anchor_domains.txt`).
* The Adjacency Index (AI): Estimates the similarity in terms of proximal domain
content by calculating the ratio between the distinct shared and distinct
unshared adjacent domains (without taking order into account)

# How does it work

BiG-SCAPE tries to (recursively) read all the GenBank files from the input 
folder (which, preferrably, correspond to identified gene clusters with a tool 
like [antiSMASH](https://antismash.secondarymetabolites.org/)). If the user has 
different subfolders in the main input directory, these can even be treated as 
different samples (and BiG-SCAPE can generate specific network files for this; 
activate with `--samples`).

BiG-SCAPE then uses the Pfam database and `hmmscan` from the HMMER (v3.1b2) 
suite to predict Pfam domains in each sequence.

For every pair of BGCs in the set(s), the pairwise distance between this BGCs is
calculated as the weighted combination of the Jaccard, AI and DDS indices. 
Network files are generated containing a number of information: the name of the 
BGCs, the raw distance between them, and data from the the three indices' 
scores. This is done taking into account different cutoff values for the 
distances (i.e. only pairs with Raw Distance < `cutoff` are written in the 
final `.network` file).

The distances for each cutoff value will be used in a clustering algorithm to
try to define 'Gene Cluster Families' (GCFs).

By default, BiG-SCAPE uses the `/product` information of antiSMASH-processed
GenBank files to separate the analysis into eight BiG-SCAPE classes: PKS Type I,
PKS Other types, NRPS, PKS/NRPS Hybrids, Saccharides, Terpenes, RiPPs and 
Others. Each have different (tuned) sets of weights for the distance components.
You can also choose to combine all BGC classes in one network file (`--mix`) and
deactivate the default classification (`--no_classify`).

See the full options with `python bigscape.py -h`.

# How to run BiG-SCAPE

## Requirements

Packages can be installed manually but using a virtual environment is
recommended. For a quick guide, see [here](Installation Guide.md)

* Python 2
* The [HMMER suite](http://hmmer.org/)
* The (processed) Pfam database. For this, download the latest `Pfam-A.hmm.gz`
file from the [Pfam website](http://pfam.xfam.org/), uncompress it and process
it using the `hmmpress` command.
* For sequence alignment (DDS score), BiG-SCAPE uses the `hmmalign` command from
the HMMER suite by default, but you can also select 
[MAFFT](http://mafft.cbrc.jp/alignment/software/) (activate with `--use_mafft`)
* Biopython
* Numpy
* scipy
* [pySAPC](https://github.com/bioinfocao/pysapc) (Affinity Propagation 
clustering algorithm with support for sparse matrices)

### Workflow

* Parses GenBank files (.gbk) and extracts CDS per BGC (fasta/*.fasta)
* Predicts domains per BGC (.domtable)
* Writes list of domains per BGC (.pfs)
* Writes selected information from filtered domtable files per BGC (.pfd)
* Writes list of sequences per domain (domains/*.fasta)
* Creates 'Arrower'-like figures for each BGC (.svg)
* Saves dictionary with list of specific domains per BGC 
(`<output dir>/BGCs.dict`)
* Calculates multiple alignments for each domain sequence file (domains/*.algn)
* Calculates distance between BGCs
* Generates network files (.network)
* Generates Network Annotation files (.tsv) with information about the input 
data
* Generates GCF labels from the clustering algorithm (.tsv)
* Generates json files for built-in visualization (.js) (*work in progress*)

BiG-SCAPE will try to re-use some of these files to continue in case the 
analysis stops (so take this into account if you e.g. change the version of the 
Pfam database or the alignment method).
