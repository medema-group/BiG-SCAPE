In order to run BiG-SCAPE, three files are needed: `bigscape.py`, `functions.py`
and `munkres.py`.

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
together to produce a compound. These proteins are most likely the most important factor when it comes to the final structure of the compound. Thus, a good distance metric should use information on the similarity of the proteins between two BGCs.

In this project, three indices are combined to define a final distance metric 
between any given pair of BGCs:
* The Jaccard index: The ratio between the distinct shared and distinct unshared
domains between two BGCs.
* The DDS (Domain Duplicate Score): Measures the similarity between two sets of 
shared domains, for each domain. If the `seqdist` distance method is selected (default), this also takes into account the sequence similarity between shared domains (calculated from multiple alignments from all sequences with the same predicted domain using MAFFT). For cases where more than one copy of a certain domain exists on at least one of the BGCs, the Hungarian algorithm is used to select the most similar Pfam domains (`munkres.py`). A special weight can be given to marked domains annotated as "anchor domains".
* The Goodman-Kruskal gamma function: Estimates the similarity of the order of 
two distinct domain sets accounting for the number of same-order and reversed pairs of domains within a certain domain count (parameter `nbhood`).


# How does it work

BiG-SCAPE tries to read all the GenBank files from the input folder (which, 
preferrably, correspond to identified gene clusters with a tool like 
[antiSMASH](https://antismash.secondarymetabolites.org/)). If the user has 
different subfolders in the main input directory, these can be treated as 
different samples (and BiG-SCAPE can generate specific network files for this; 
see the parameters with the `--help` option).

BiG-SCAPE then uses `hmmscan` from the HMMER (v3.1b2) package to predict Pfam 
domains in each sequence.

For every pair of BGCs in the set(s), the pairwise distance between this BGCs is
calculated as the weighted combination of the Jaccard, GK and DDS indices. Network files are generated containing a number of information: the name of the BGCs, the raw distance, the three indices' score, a -log2 score of the raw distance and more. This is done taking into account different cutoff values for the distances (i.e. only rows with raw distance <= cutoff are written in the file).


# How to run/re-run BiG-SCAPE

Before using BiG-SCAPE, make sure that `MAFFT` and `hmmscan` are available in 
your path. `hmmscan` also uses the Pfam library (specifically, Pfam-A.hmm formatted with the `hmmpress` command).

Run as `python bigscape.py <options>`.

There are some options to re-run BiG-SCAPE and skip some of the most 
time-consuming calculations:

## Parameter (nothing)

### What it does
* Parses GenBank files (.gbk) and extracts CDS per BGC (.fasta)
* Predicts domains per BGC (.domtable)
* Writes list of sequences per domain (.fasta)
* Writes list of domains per BGC (.pfs)
* Writes selected information from domtable files per BGC (.pfd)
* Saves dictionary with list of specific domains per BGC 
(`<output dir>/BGCs.dict`)
* When using the `seqdist` distance method, calculates multiple alignments using
MAFFT (.algn) for domains represented by more than one sequence. Saves alignment information in `<output dir>/DMS.dict`
* Calculates distance between BGCs
* Generates network files

### When to use it
* First time
* When changing fundamental parameters like `domain_overlap_cutoff` or 
`min_bgc_size`. Note that BiG-SCAPE only works with files which do not yet have their corresponding processed output (so it's possible to save time when resuming from an unexpected end of the run). This means that re-running BiG-SCAPE with different values for these parameters will not have an effect unless you empty your output directory or choose a new one.
* When you are changing versions of the pfam database (if you have previous results from past runs, you can use the `--force_hmmscan` parameter)

### What you need
* A list of GenBank files (.gbk). If either the `seqdist` or `domain_dist` 
methods are selected with the "S" option, BiG-SCAPE will treat subfolders within
your input directory as "samples". The name of each sample is taken from each 
subfolder with .gbk files and, as it is used in the final network files duplicated subfolder/sample names is not allowed.


## Parameter `--skip_mafft`

### What it does
Skips domain sequence alignment using `MAFFT`. Oherwise:
* Calculates distance between BGCs
* Generates network files

### When to use it
* You already ran BiG-SCAPE once
* You are using the `seqdist` distance method
* You want to change the `nbhood` parameter for the Goodman-Kruskal index
* When you had only previously calculated distance within samples (i.e. not 
all-vs-all distances)
* When you want to change the distance method from first run (e.g. from 
`domaindist` to `seqdist`)
* When you need to recalculate distances (e.g. to change the anchor domains 
list)

### What you need
* The original .gbk files. These are necessary to rebuild the original structure
* The list of domains per BGC (.pfs files). These are necessary for distance 
calculations
* The output from `hmmscan` (.domtable files). These are necesary to see whether
a particular .gbk file had no predicted domains (in which case it'll be taken out of the analysis)
* The BGCs dictionary (`<output dir>/BGCs.dict`). This contains a list of 
domains per BGC (and for each of these domains, a list of all predicted locations within the BGC). Needed for distance calculations
* The DMS dictionary (`<output dir>/DMS.dict`). This contains information from 
the sequence alignment

## Parameter `--skip_all`

### What it does
* Skips all calculations
* Generates new network files

### When to use it
* When you only want to generate new network files: use different cutoffs, 
include/exclude singletons, use input file structure to generate network files 
per sample or change weights for the indices

### What you need
* The output from `hmmscan`.
* The network file with the appropriate distance method and complete set of 
distances. This means either `networkfile_seqdist_all_vs_all_c1.network` for 
the `seqdist` method or `networkfile_domain_dist_all_vs_all_c1.network` for 
`domain_dist`. BiG-SCAPE will make use of these pre-calculated distances so 
these cutoff=1 value is added by default in the list of cutoffs
