import glob
import logging
import os
import pyhmmer
from src.data.bigslice import get_antismash_domains, get_bigslice_subset

from src.data import Database


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
            if accession in bigslice_accessions:
                model_type = 1
                name = "AS-" + name
                bigslice_accessions.remove(accession)
            
            if name in bigslice_accessions:
                model_type = 1
                name = "AS-" + name
                bigslice_accessions.remove(name)
                

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
    bigslice_pfam_path = os.path.join(
        run.bigslice.bigslice_data_path,
        "biosynthetic_pfams",
        "Pfam-A.biosynthetic.hmm"
    )

    with pyhmmer.plan7.HMMFile(bigslice_pfam_path) as hmm_file:
        for profile in hmm_file:
            profile: pyhmmer.plan7.Profile

            # skip anything that doesn't have an accession
            if profile.accession is None:
                continue

            accession = profile.accession.decode()
            # only check ones prefixed with AS-
            if not accession.startswith("AS-"):
                continue

            # then check if it's there. we need to trim the "AS-" for this
            if accession[3:] not in antismash_domains:
                continue

            model_type = 2
            accession = None
            if profile.accession is not None:
                accession = profile.accession.decode()
            name = profile.name.decode()

            hmm_id = database.insert("hmm", {
                "accession": accession,
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
                name = profile.name.decode()

                if profile.accession is not None:
                    accession = profile.accession.decode()
                else:
                    accession = name



                hmm_id = database.insert("hmm", {
                    "accession": accession,
                    "name": name,
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
