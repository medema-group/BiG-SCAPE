from os import path, makedirs, remove, rename, SEEK_END, sched_getaffinity
from shutil import copy, rmtree, copyfileobj
from hashlib import md5
import urllib.request
import gzip
import csv
import glob
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
            return(1)

        print("Extracting bigslice_models.tar.gz...")

        with tarfile.open(zipped_file, "r:gz") as fp:
            fp.extractall(path=models_folder)

        print("done! (please remove the downloaded tar.gz file manually)")
        return(0)

    else:
        print("models folder exists!")
        return(1)


def md5sum(filename):
    hash = md5()
    with open(filename, "rb") as f:
        for chunk in iter(lambda: f.read(128 * hash.block_size), b""):
            hash.update(chunk)
    return hash.hexdigest()
