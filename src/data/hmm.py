import glob
import logging
import os
import pyhmmer
from src.big_scape.run.base import Run

from src.data import Database

def get_bigslice_subset():
    """Returns a list of pfam accessions that are included in the bigslice subset of
    accessions, and which will be used later on to do a pre-filter step
    """
    bigslice_accessions = set()
    bigslice_names = set()

    biopfam_path = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
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
    as_domain_info_file = os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
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


def load_hmms(run, database: Database):
    """Loads the hmm information into the database"""

    # main pfam
    hmm_file_path = os.path.join(run.directories.pfam, "Pfam-A.hmm")
    logging.info(" Loading HMM information into database")

    # get list of accessions and descriptions which bigslice uses from raw Pfam_A
    bigslice_accessions, bigslice_descriptions = get_bigslice_subset()
    
    # get list of antismash domains
    antismash_domains, antismash_domain_info = get_antismash_domains(run)

    # later we need to know which subpfams belong to which domains
    # in order to do this we need to keep track of all accession <-> id relations
    domain_accession_ids = {}
    # sometimes the accession is unavailable and we need to use the names instead
    domain_name_ids = {}

    # add all pfam domains
    with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
        for profile in hmm_file:
            profile: pyhmmer.plan7.Profile

            model_type = 0
            accession = profile.accession.decode()
            name = profile.name.decode()

            # if bigslice prefiltering is enabled we need to check whether
            # domains are parent antismash domains that are later used
            # to calculate features for subpfam domains. to do this we
            # prefix names with AS- if it is a parent antismash domain
            # needed for bigslice feature calculation and not a problem for
            # BiG-SCAPE
            if accession in bigslice_accessions or name in bigslice_accessions:
                model_type = 1
                name = "AS-" + name

            hmm_id = database.insert("hmm", {
                "accession": accession,
                "name": name,
                "model_length": profile.M,
                "model_type": model_type
            }, True)

            domain_accession_ids[accession] = hmm_id
            domain_name_ids[name] = hmm_id
    database.commit_inserts()

    # end early if bigslice prefilter is disabled
    if not run.bigslice.use_bigslice:
        return

    # add any missing antismash domains
    # base data dir where antismash files were extracted to
    antismash_data_dir = os.path.join(
    run.bigslice.bigslice_data_path,
    "antismash",
    run.bigslice.antismash_tar_folder,
    "antismash",
    "detection",
    "hmm_detection",
    "data")

    for antismash_domain in antismash_domains:
        # get domain info
        domain_info = antismash_domain_info[antismash_domain]

        hmm_file_path = os.path.join(antismash_data_dir, domain_info["filename"])
        with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
            # this should only ever be one profile, but just to be sure let's loop through
            for profile in hmm_file:
                profile: pyhmmer.plan7.Profile

                model_type = 2
                accession = None
                if profile.accession is not None:
                    accession = profile.accession.decode()
                name = "AS-" + profile.name.decode()

                hmm_id = database.insert("hmm", {
                    "accession": name,
                    "name": name,
                    "model_length": profile.M,
                    "model_type": model_type
                }, True)

                domain_name_ids[name] = hmm_id
    database.commit_inserts()

    # finally, sub pfams
    logging.info("Loading BiG-SLICE subpfams into database")
    sub_pfam_path = os.path.join(
        run.bigslice.bigslice_data_path,
        "sub_pfams",
        "hmm"
    )
    for hmm_path in glob.iglob(os.path.join(
        sub_pfam_path,
        "*.hmm")):
        with pyhmmer.plan7.HMMFile(hmm_path) as hmm_file:
            for profile in hmm_file:
                profile: pyhmmer.plan7.Profile

                accession = None
                if profile.accession is not None:
                    accession = profile.accession.decode()


                hmm_id = database.insert("hmm", {
                    "accession": accession,
                    "name": profile.name.decode(),
                    "model_length": profile.M,
                    "model_type": 3
                }, True)

                # add relations
                base_name = hmm_path.rsplit("/", maxsplit=1)[1].split(".subpfams.hmm")[0]
                
                parent_id = None
                # first check if a HMM with a corresponding name is found
                if base_name in domain_name_ids:
                    parent_id = domain_name_ids[base_name]
                
                # we can check whether a version with

                # if not, then we check if there is an accession
                elif base_name in domain_accession_ids:
                    parent_id = domain_accession_ids[base_name]

                # otherwise we will just have to skip it
                else:
                    continue

                database.insert("subpfam", {
                    "hmm_id": hmm_id,
                    "parent_hmm_id": parent_id
                })
        database.commit_inserts()

def from_id(database: Database, hmm_id: int):
    """Returns a hmm row based on a given ID"""
    return database.select("hmm",
                           f"where id = {hmm_id}",
                           props="*")

def from_accession(database: Database, accession: str):
    """Returns a hmm row based on a given accession"""
    rows = database.select("hmm",
                           f"where accession = \"{accession}\"",
                           props="*")
    if len(rows) == 0:
        return None
    return rows[0]

def from_model_type(database: Database, model_type: int):
    """Returns a list of hmm details from a given model type
    Inputs:
        - model types:
            0: only found in BiG-SCAPE analysis
            1: found in BiG-SCAPE and BiG-SLICE analysis
            2: only found in BiG-SLICE
    """
    rows = database.select(
        "hmm",
        f"where model_type = \"{model_type}\""
    )
    return rows
