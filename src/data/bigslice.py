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


def generate_core_bio_hmm(run):
    """temp"""
    dir_path = run.bigslice.bigslice_data_path
    tmp_dir_path = path.join(dir_path, "tmp")

    # create temporary directory
    if not path.exists(tmp_dir_path):
        makedirs(tmp_dir_path)

    # prepare generate_hmm steps

    antismash_folder = path.join(tmp_dir_path,
                                 run.bigslice.ANTISMASH_VERSION,
                                 "antismash")

    biosyn_pfam_tsv = path.join(dir_path, "biosynthetic_pfams", "biopfam.tsv")
    biosyn_pfam_hmm = path.join(
        dir_path, "biosynthetic_pfams", "Pfam-A.biosynthetic.hmm")

    sub_pfams_tsv = path.join(dir_path, "sub_pfams", "corepfam.tsv")
    sub_pfams_hmms = path.join(dir_path, "sub_pfams", "hmm")

    def get_pressed_hmm_filepaths(hmm_file_path):
        base_name = hmm_file_path
        return (
            base_name + ".h3i",
            base_name + ".h3f",
            base_name + ".h3m",
            base_name + ".h3p",
        )

    # create subpfam directory
    if not path.exists(sub_pfams_hmms):
        makedirs(sub_pfams_hmms)

    # download and extract antiSMASH models
    if not path.exists(antismash_folder):
        antismash_zipped_file = path.join(tmp_dir_path, "antismash.tar.gz")
        if not path.exists(antismash_zipped_file):
            print("Downloading antismash.tar.gz...")
            urllib.request.urlretrieve(
                run.bigslice.ANTISMASH_URL, path.join(
                    tmp_dir_path, "antismash.tar.gz"))

        print("Extracting antismash.tar.gz...")
        with tarfile.open(antismash_zipped_file, "r:gz") as as_zipped:
            as_zipped.extractall(path=tmp_dir_path)

    # parse antiSMASH hmmdetails.txt
    antismash_domains = {}
    antismash_from_pfams = set()
    antismash_data_folder = path.join(antismash_folder,
                                      "detection",
                                      "hmm_detection",
                                      "data")
    hmmdetails_path = path.join(antismash_data_folder,
                                "hmmdetails.txt")
    with open(hmmdetails_path, "r") as hmmdetails:
        for line in hmmdetails:
            line = line.rstrip()
            name, desc, cutoff, filename = line.split("\t")
            cutoff = int(cutoff)
            acc = None

            with open(path.join(
                    antismash_data_folder,
                    filename), "r") as ashmm:
                for line in ashmm:
                    line = line.rstrip()
                    if line.startswith("ACC "):
                        acc = line.split(" ")[-1]
                    elif line.startswith("HMM "):
                        break

            antismash_domains[name] = {
                "desc": desc,
                "cutoff": cutoff,
                "filename": filename
            }

            if acc and acc.startswith("PF"):
                antismash_from_pfams.add(acc.split(".")[0])

    # parse antiSMASH cluster_rules to get a list of
    # core domains
    cluster_rules_dir = path.join(antismash_folder,
                                  "detection",
                                  "hmm_detection",
                                  "cluster_rules")
    antismash_core_list = set()
    for cluster_rules_file in glob.iglob(
            path.join(cluster_rules_dir, "*.txt")):
        with open(cluster_rules_file, "r") as cr_handle:
            rules = parse_antismash_rules(cr_handle)
            for rule, rule_prop in rules.items():
                antismash_core_list.update(
                    fetch_antismash_domain_names(rule_prop["conditions"]))

    as_domains = set(antismash_domains.keys())
    for core in antismash_core_list:
        if core not in as_domains:
            raise Exception(core + " is not found in antiSMASH domains!")


    # (down)loads Pfam-A.hmm
    if not path.exists(path.join(tmp_dir_path, "Pfam-A.hmm.gz")):
        print("Downloading Pfam-A.hmm.gz...")
        urllib.request.urlretrieve(
            run.bigslice.PFAM_DB_URL, path.join(
                tmp_dir_path, "Pfam-A.hmm.gz"))

    # load biosynthetic pfams list
    biosynthetic_pfams = []
    with open(biosyn_pfam_tsv, "r") as biopfam_tsv:
        reader = csv.DictReader(biopfam_tsv, dialect="excel-tab")
        for row in reader:
            if row["Status"] == "included":
                # check if included in antismash
                if row["Acc"].split(".")[0] in antismash_from_pfams:
                    print(row["Acc"] + " is in antiSMASH pfam, skipping..")
                    continue
                biosynthetic_pfams.append(row["Acc"])

    # apply biosynthetic pfams filtering
    with gzip.open(path.join(tmp_dir_path, "Pfam-A.hmm.gz"), "rt") as pfam:
        with open(biosyn_pfam_hmm, "w") as biopfam:
            print("Generating Pfam-A.biosynthetic.hmm...")
            temp_buffer = ""  # for saving a temporary hmm entry
            skipping = False
            for line in pfam:
                if line.startswith("//") and len(temp_buffer) > 0:  # flush
                    if not skipping:
                        biopfam.write(temp_buffer)
                    temp_buffer = ""
                    skipping = False
                if skipping:
                    continue
                temp_buffer += line
                if line.startswith("ACC "):
                    pfam_acc = line.split(" ")[-1].rstrip()
                    try:
                        biosynthetic_pfams.remove(pfam_acc)
                    except ValueError:
                        skipping = True
            biopfam.write("//\n")

        assert len(biosynthetic_pfams) == 0

        # add antismash domains
        with open(biosyn_pfam_hmm, "a") as biopfam:
            for as_name, as_data in antismash_domains.items():
                with open(path.join(
                        antismash_data_folder,
                        as_data["filename"]), "r") as ashmm:
                    for line in ashmm:
                        line = line.rstrip()
                        if line.startswith("NAME "):
                            line = "NAME  AS-" + as_name + "\n"
                            line += "ACC   AS-" + as_name
                        elif line.startswith("ACC "):
                            continue
                        elif line.startswith("DESC "):
                            line = "DESC  " + as_data["desc"]
                        elif line.startswith(("TC ", "GA ", "NC ")):
                            continue
                        elif line.startswith("STATS LOCAL MSV "):
                            cutoff_str = str(as_data["cutoff"]) + ".0"
                            biopfam.write("GA    {} {};\n".format(
                                cutoff_str, cutoff_str))
                            biopfam.write("TC    {} {};\n".format(
                                cutoff_str, cutoff_str))
                            biopfam.write("NC    {} {};\n".format(
                                cutoff_str, cutoff_str))
                        biopfam.write(line + "\n")

    hmm_presses = get_pressed_hmm_filepaths(biosyn_pfam_hmm)
    hmm_pressed = True
    for hmm_press_file in hmm_presses:
        if not path.exists(hmm_press_file):
            hmm_pressed = False
            break
    if not hmm_pressed:
        print("Running hmmpress on Pfam-A.biosynthetic.hmm")
        for hmm_press_file in hmm_presses:
            if path.exists(hmm_press_file):
                remove(hmm_press_file)
        if subprocess.call([
            "hmmpress",
            biosyn_pfam_hmm
        ]) > 0:
            return

    # ----------- sub_pfams ----------- #

    # fetch subpfams list
    sub_pfams = {}
    with open(sub_pfams_tsv, "r") as corepfam:
        corepfam.readline()  # skip header
        for line in corepfam:
            [pfam_accession, pfam_name, pfam_desc] = line.rstrip().split("\t")
            if pfam_accession in sub_pfams:
                raise Exception(
                    "Duplicated core pfam found in " + sub_pfams_tsv)
            else:
                # check if included in antismash
                if pfam_accession.split(".")[0] in antismash_from_pfams:
                    print(pfam_accession +
                          " is in antiSMASH pfam, skipping subpfam..")
                    continue
                sub_pfams[pfam_accession] = {
                    "name": pfam_name,
                    "desc": pfam_desc
                }
    for as_core in antismash_core_list:  # from antiSMASH
        if as_core in antismash_domains:
            acc = "AS-" + as_core
            if acc in sub_pfams:
                raise Exception(
                    "Duplicated core pfam found: " + acc)
            sub_pfams[acc] = {
                "name": "AS-" + as_core,
                "desc": antismash_domains[as_core]["desc"]
            }
        else:
            raise Exception(as_core + " is not in antiSMASH domain list!")

    # generate core pfams-only HMM
    core_hmms_path = path.join(tmp_dir_path, "core_pfams.hmm")
    with open(biosyn_pfam_hmm, "r") as full_hmm:
        with open(core_hmms_path, "w") as core_hmm:
            buf_text = ""
            cur_acc = ""
            cur_leng = 0
            have_ga = False
            for line in full_hmm:
                buf_text += line
                if line.startswith("ACC"):
                    cur_acc = line.rstrip().split(" ")[-1]
                elif line.startswith("GA"):
                    have_ga = True
                elif line.startswith("LENG"):
                    cur_leng = int(
                        line.rstrip().split(" ")[-1])
                elif line.startswith("//"):
                    if cur_acc in sub_pfams:
                        sub_pfams[cur_acc]["leng"] = cur_leng
                        if not have_ga:
                            before, after = buf_text.split("CKSUM")
                            buf_text = before + \
                                "GA    20.00 20.00;\n" + \
                                "CKSUM" + after
                        core_hmm.write(buf_text)
                    buf_text = ""
                    cur_acc = ""
                    have_ga = False
            if cur_acc in sub_pfams:
                if not have_ga:
                    before, after = buf_text.split("CKSUM")
                    buf_text = before + \
                        "GA    20.00 20.00;\n" + \
                        "CKSUM" + after
                core_hmm.write(buf_text)

    # download reference proteins dataset
    ref_prot_filename = path.basename(run.bigslice.REFERENCE_PROTEINS_URL)
    stored_ref_prot_filename = "subpfam_refprot.fa"
    if not path.exists(path.join(tmp_dir_path, stored_ref_prot_filename)):
        if not path.exists(path.join(tmp_dir_path, ref_prot_filename)):
            print("Downloading " + ref_prot_filename)
            urllib.request.urlretrieve(
                run.bigsliceREFERENCE_PROTEINS_URL, path.join(
                    tmp_dir_path, ref_prot_filename))
        file_ext = path.splitext(ref_prot_filename)[1]
        if file_ext in ["fasta", "fa"]:
            rename(
                path.join(tmp_dir_path, ref_prot_filename),
                path.join(tmp_dir_path, stored_ref_prot_filename)
            )
        elif file_ext == ".gz":
            print("Extracting reference proteins...")
            with gzip.open(
                path.join(tmp_dir_path, ref_prot_filename), 'rb'
            ) as f_in:
                with open(
                    path.join(tmp_dir_path, stored_ref_prot_filename), 'wb'
                ) as f_out:
                    copyfileobj(f_in, f_out)
        else:
            raise Exception("Unrecognized file format! " + file_ext)

    print("HMMScanning reference proteins...")

    # run hmmscan to get aligned fasta files
    ref_prot_hmmtxt = path.join(
        tmp_dir_path, "subpfam_refprot_hits.txt")
    if not path.exists(ref_prot_hmmtxt):
        command = [
            "hmmsearch",
            "--acc",
            "--cut_ga",
            "--cpu", str(len(sched_getaffinity(0))),
            "-o", ref_prot_hmmtxt,
            core_hmms_path,
            path.join(tmp_dir_path, stored_ref_prot_filename)
        ]
        ret = subprocess.run(command, check=True)
        if ret.returncode != 0:
            raise Exception("Error doing hmmsearch")
    else:
        with open(ref_prot_hmmtxt, "r") as fp:
            if next(reversed_fp_iter(fp)).rstrip() != "[ok]":
                # file is broken, remove
                raise Exception(
                    ref_prot_hmmtxt +
                    " is broken, please remove and rerun the script")

    print("Parsing hmmscan results")
    # parse hmmtxt into alignment fastas
    core_pfam_hit_counts = {acc: 0 for acc in sub_pfams.keys()}
    for run_result in parse(ref_prot_hmmtxt, "hmmer3-text"):
        pfam_acc = run_result.accession
        model_len = sub_pfams[pfam_acc]["leng"]
        aligned_multifasta_path = path.join(
            tmp_dir_path,
            pfam_acc + ".aligned.fa"
        )

        with open(aligned_multifasta_path, "w") as fa:
            for hsp in run_result.hsps:
                # get hits aligned to model
                aligned_to_model = ""
                query_seq = str(hsp.query.seq)
                hit_seq = str(hsp.hit.seq)
                for i in range(hsp.query_start):
                    aligned_to_model += "-"
                for i in range(len(query_seq)):
                    if query_seq[i] != ".":
                        aligned_to_model += hit_seq[i]
                for i in range(model_len - len(aligned_to_model)):
                    aligned_to_model += "-"

                # write to fasta
                fa_acc = hsp.hit.id + "|" + \
                    str(hsp.hit_start) + "-" + str(hsp.hit_end)
                fa.write(">{}\n{}\n".format(fa_acc, aligned_to_model))
                core_pfam_hit_counts[pfam_acc] += 1

    # filter for subpfams not having any hits
    for acc, count in core_pfam_hit_counts.items():
        if count < run.bigslice.MIN_PROTEIN_SEQUENCES:
            print("skipping subpfam " + acc + ", not enough reference...")
            del sub_pfams[acc]

    # build subpfams
    for pfam_accession, pfam_properties in sub_pfams.items():
        subpfam_hmm_path = path.join(
            sub_pfams_hmms, "{}.subpfams.hmm".format(pfam_accession))
        if not path.exists(subpfam_hmm_path):
            print("Building {}...".format(subpfam_hmm_path))
            aligned_multifasta_path = path.join(
                tmp_dir_path,
                pfam_accession + ".aligned.fa"
            )
            temp_hmm_path = path.splitext(aligned_multifasta_path)[
                0] + ".subpfams.hmm"
            if not path.exists(temp_hmm_path):
                build_subpfam(
                    run,
                    aligned_multifasta_path,
                    temp_hmm_path
                )
            copy(temp_hmm_path, subpfam_hmm_path)

        # check hmmpress
        hmm_presses = get_pressed_hmm_filepaths(subpfam_hmm_path)
        hmm_pressed = True
        for hmm_press_file in hmm_presses:
            if not path.exists(hmm_press_file):
                hmm_pressed = False
                break
        if not hmm_pressed:
            print("Running hmmpress on {}".format(subpfam_hmm_path))
            for hmm_press_file in hmm_presses:
                if path.exists(hmm_press_file):
                    remove(hmm_press_file)
            if subprocess.call([
                "hmmpress",
                subpfam_hmm_path
            ]) > 0:
                return


def build_subpfam(run, input_fasta, output_hmm):
    """Builds subpfam hmm file"""
    # fetch sequences into numpy array
    labels = []
    sequences = []
    with open(input_fasta, "r") as fa:
        for line in fa:
            line = line.rstrip()
            if len(line) > 0:
                if line.startswith(">"):
                    labels.append(line[1:])
                else:
                    if line.count("-") > (len(line) * 0.75):
                        if len(labels) > 0:
                            del labels[-1]
                    else:
                        sequences.append(list(line))
    assert len(labels) == len(sequences)
    uniques = np.unique(sequences)
    X = np.searchsorted(uniques, sequences)

    # calculate pdist
    dist = pairwise_distances(
        X, metric='hamming', n_jobs=-1)

    # clustering
    cl = AgglomerativeClustering(
        n_clusters=min(run.bigslice.MAX_CLADE_NUMBERS, max(
            run.bigslice.MIN_CLADE_NUMBERS, int(X.shape[0] / 200))),
        affinity='precomputed',
        linkage='complete')
    cl.fit(dist)

    clades = {}
    for i, label in enumerate(cl.labels_):
        if label not in clades:
            clades[label] = []
        clades[label].append(i)

    # hmmbuild
    with TemporaryDirectory() as temp_dir:
        pfam_name = path.splitext(path.basename(input_fasta))[0]
        for cl, clade in clades.items():
            hmm_name = "{}_c{}".format(pfam_name, cl)
            with open(path.join(
                    temp_dir, "{}.fa".format(hmm_name)), "w") as fa:
                for i in clade:
                    seq = "".join([c for c in sequences[i]])
                    fa.write(">{}\n{}\n".format(labels[i], seq))
            command = [
                "hmmbuild",
                "--cpu",
                str(len(sched_getaffinity(0))),
                "-n",
                hmm_name,
                "-o",
                "/dev/null",
                path.join(temp_dir, "{}.hmm".format(hmm_name)),
                path.join(temp_dir, "{}.fa".format(hmm_name))
            ]
            subprocess.call(" ".join(command), shell=True)
        with open(output_hmm, "w") as hm:
            for cl in clades:
                hmm_name = "{}_c{}".format(pfam_name, cl)
                with open(path.join(
                        temp_dir, "{}.hmm".format(hmm_name)), "r") as sm:
                    for line in sm.readlines():
                        hm.write(line)

def parse_antismash_rules(file_handle):
    """Parses antismash rules file"""
    conds = {}
    cur_cond = {}
    parsing_cond = False
    for line in file_handle:
        line = line.split("#")[0]
        line = line.rstrip().lstrip()
        if line.startswith("RULE "):
            if "rule" in cur_cond:
                conds[cur_cond["rule"]] = cur_cond
            cur_cond = {
                "rule": line.split(" ")[-1]
            }
            parsing_cond = False
        elif parsing_cond:
            if len(line) > 0:
                cur_cond["conditions"] += " " + line
        elif line.startswith("COMMENT "):
            cur_cond["comment"] = line.split(" ")[-1]
        elif line.startswith("CUTOFF "):
            cur_cond["cutoff"] = int(line.split(" ")[-1])
        elif line.startswith("NEIGHBOURHOOD "):
            cur_cond["neighbourhood"] = int(line.split(" ")[-1])
        elif line.startswith("CONDITIONS "):
            cur_cond["conditions"] = line.split("CONDITIONS ")[-1]
            parsing_cond = True
    if "rule" in cur_cond:
        conds[cur_cond["rule"]] = cur_cond

    return conds


def fetch_antismash_domain_names(conditions):
    """Fetches antismash domain names"""
    # replace all non-token with space
    non_tokens = [
        "minimum(", "minscore(",
        " and ", "not ",
        " or ", "cds(",
        "[", "]", ",", ")", "("  # shorter chars should be in the back
    ]
    for nt in non_tokens:
        conditions = conditions.replace(nt, " ")
    return set([
        token for token in conditions.split(" ")
        if len(token) > 0 and not token.isdigit()
    ])


def reversed_fp_iter(fp, buf_size=8192):
    """a generator that returns the lines of a file in reverse order
    ref: https://stackoverflow.com/a/23646049/8776239
    """
    segment = None
    offset = 0
    fp.seek(0, SEEK_END)
    file_size = remaining_size = fp.tell()
    while remaining_size > 0:
        offset = min(file_size, offset + buf_size)
        fp.seek(file_size - offset)
        buffer = fp.read(min(remaining_size, buf_size))
        remaining_size -= buf_size
        lines = buffer.splitlines(True)
        # the first line of the buffer is probably not a complete line so
        # we'll save it and append it to the last line of the next buffer
        # we read
        if segment is not None:
            # if the previous chunk starts right from the beginning of line
            # do not concat the segment to the last line of new chunk
            # instead, yield the segment first
            if buffer[-1] == '\n':
                # print 'buffer ends with newline'
                yield segment
            else:
                lines[-1] += segment
        segment = lines[0]
        for index in range(len(lines) - 1, 0, -1):
            if len(lines[index]):
                yield lines[index]
    # Don't yield None if the file was empty
    if segment is not None:
        yield segment

