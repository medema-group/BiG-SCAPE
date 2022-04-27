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
