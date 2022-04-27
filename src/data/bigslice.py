import glob
from os import path, makedirs, remove, rename, SEEK_END, sched_getaffinity
from shutil import copy, rmtree, copyfileobj
from hashlib import md5
import urllib.request
import gzip
import csv
from tempfile import TemporaryDirectory
import subprocess
import tarfile

# external imports
import numpy as np
import pyhmmer
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import pairwise_distances
from Bio.SearchIO import parse


def download_bigslice_db(run):
    """Downloads the necessary files for bigslice feature generation"""
    models_folder = run.bigslice.bigslice_data_path

    # if not path.exists(models_folder):
    if not path.exists(path.join(models_folder, "biosynthetic_pfams")):
        zipped_file = "bigslice_models.tar.gz"
        if not path.exists(zipped_file):
            print("Downloading bigslice_models.tar.gz...")
            urllib.request.urlretrieve(run.bigslice_db_location, zipped_file)

        print("Checking MD5 sums...")
        md5sum_downloaded = md5sum(zipped_file)
        if run.bigslice.bigslice_db_md5 != md5sum_downloaded:
            print("'{}' vs '{}': doesn't match!".format(
                run.bigslice.bigslice_db_md5, md5sum_downloaded))

        print("Extracting bigslice_models.tar.gz...")

        with tarfile.open(zipped_file, "r:gz") as fp:
            fp.extractall(path=models_folder)

    else:
        print("models folder exists!")


def download_antismash_files(run):
    """download and extract antiSMASH models"""
    antismash_folder = path.join(
        run.bigslice.bigslice_data_path,
        "antismash"
    )
    if not path.exists(antismash_folder):
        antismash_zipped_file = path.join(
            run.bigslice.bigslice_data_path,
            "antismash.tar.gz"
        )
        if not path.exists(antismash_zipped_file):
            print("Downloading antismash.tar.gz...")
            urllib.request.urlretrieve(
                run.bigslice.antismash_url, antismash_zipped_file)

        print("Extracting antismash.tar.gz...")
        with tarfile.open(antismash_zipped_file, "r:gz") as as_zipped:
            as_zipped.extractall(path=antismash_folder)
        print("Done extracting antismash.tar.gz")


def md5sum(filename):
    hash = md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(128 * hash.block_size), b""):
            hash.update(chunk)
    return hash.hexdigest()


def get_bigslice_subset():
    """Returns a list of pfam accessions that are included in the bigslice subset of
    accessions, and which will be used later on to do a pre-filter step
    """
    bigslice_accessions = set()
    bigslice_names = set()

    biopfam_path = path.join(
                    path.dirname(path.abspath(__file__)),
                    "biopfam.tsv"
    )

    # corepfam_path = os.path.join(
    #                 os.path.dirname(os.path.abspath(__file__)),
    #                 "corepfam.tsv"
    # )

    with open(biopfam_path, encoding="utf-8") as bio_pfam_file:
        for line in bio_pfam_file:
            lineparts = line.rstrip().split("\t")

            if lineparts[3] == "included":
                bigslice_accessions.add(lineparts[0])
                bigslice_names.add(lineparts[1])

    # with open(corepfam_path, encoding="utf-8") as core_pfam_file:
    #     for line in core_pfam_file:
    #         lineparts = line.rstrip().split("\t")

    #         bigslice_accessions.add(lineparts[0])
    #         bigslice_names.add(lineparts[1])

    print(len(bigslice_accessions))

    return bigslice_accessions, bigslice_names

def get_antismash_domains(run):
    """Returns a set of antismash domain names and a dictionary of domain info
    """
    as_domains = set()
    as_domain_info = dict()
    as_domain_info_file = path.join(
                    path.dirname(path.abspath(__file__)),
                    "hmmdetails.txt"
    )
    for line in open(as_domain_info_file, "r"):
        name, desc, cutoff, filename = line.rstrip().split("\t")
        as_domains.add(name)
        as_domain_info[name] = {
            "desc": desc,
            "cutoff": cutoff,
            "filename": filename
        }
    return as_domains, as_domain_info


def get_bigslice_biosynth_profiles(run, bigslice_accessions):
    """Gets a list of pyhmmer profiles specific for biosynthetic hmm search
    """
    bigslice_profiles = list()

    # start with the antismash profiles that weren't included
    # in the full Pfam-A hmm

    # we'll get them from the BiG-SLICE models
    bigslice_pfam_path = path.join(
        run.bigslice.bigslice_data_path,
        "biosynthetic_pfams",
        "Pfam-A.biosynthetic.hmm"
    )

    with pyhmmer.plan7.HMMFile(bigslice_pfam_path) as hmm_file:
        optimized_profiles = hmm_file.optimized_profiles()
        # this should only ever be one profile, but just to be sure let's loop through
        for profile in optimized_profiles:
            profile: pyhmmer.plan7.OptimizedProfile

            # cheat to set the trusted cutoffs to the gathering cutoffs
            # used in BiG-SLICE
            profile.cutoffs.trusted = profile.cutoffs.gathering

            # only include ones we didn't previously include
            if profile.accession.decode() in bigslice_accessions:
                bigslice_profiles.append(profile)
    return bigslice_profiles


def get_bigslice_subpfam_profiles(run):
    """Gets a list of pyhmmer profiles for subpfam hmmsearch
    """
    bigslice_profiles = list()
    sub_pfam_path = path.join(
        run.bigslice.bigslice_data_path,
        "sub_pfams",
        "hmm"
    )
    for hmm_path in glob.iglob(path.join(
        sub_pfam_path,
        "*.hmm")):
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            optimized_profiles = hmm_file.optimized_profiles()
            for profile in optimized_profiles:
                profile: pyhmmer.plan7.OptimizedProfile

                accession = None
                name = profile.name.decode()

                if profile.accession is not None:
                    accession = profile.accession.decode()
                else:
                    accession = name

                # we can set the cutoff manually for this as it is in BiG-SLICE
                cutoff = 20
                profile.cutoffs.trusted = [cutoff, cutoff]
                
                # profile.accession = accession.encode()

                bigslice_profiles.append(profile)

    return bigslice_profiles

