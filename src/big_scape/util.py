import os
import json
import sys

from src.bgctools import sort_bgc
from src.js import add_to_bigscape_results_js
from src.pfam import get_domain_list

def get_ordered_domain_list(RUN, CLUSTER_BASE_NAMES):
    print(" Reading the ordered list of domains from the pfs files")

    DOMAIN_LIST = {} # Key: BGC. Item: ordered list of domains
    for outputbase in CLUSTER_BASE_NAMES:
        pfsfile = os.path.join(RUN.directories.pfs, outputbase + ".pfs")
        if os.path.isfile(pfsfile):
            DOMAIN_LIST[outputbase] = get_domain_list(pfsfile)
        else:
            sys.exit(" Error: could not open " + outputbase + ".pfs")
    return DOMAIN_LIST

def fetch_genome_list(run, input_clusters_idx, cluster_names, mibig_set, bgc_info, gen_bank_dict):
    genomes = []
    classes = []
    cluster_names_to_genomes = {}
    cluster_names_to_classes = {}
    for idx, bgc in enumerate(cluster_names):
        if bgc in mibig_set:
            continue
        input_clusters_idx.append(idx)
        # get class info
        product = bgc_info[bgc].product
        predicted_class = sort_bgc(product)
        if predicted_class not in classes:
            cluster_names_to_classes[bgc] = len(classes)
            classes.append(predicted_class)
        else:
            cluster_names_to_classes[bgc] = classes.index(predicted_class)
        # get identifier info
        identifier = ""
        if len(bgc_info[bgc].organism) > 1:
            identifier = bgc_info[bgc].organism
        else: # use original genome file name (i.e. exclude "..clusterXXX from antiSMASH run")
            file_name_base = os.path.splitext(os.path.basename(gen_bank_dict[bgc][0]))[0]
            identifier = file_name_base.rsplit(".cluster", 1)[0].rsplit(".region", 1)[0]
        if len(identifier) < 1:
            identifier = "Unknown Genome {}".format(len(genomes))
        if identifier not in genomes:
            cluster_names_to_genomes[bgc] = len(genomes)
            genomes.append(identifier)
        else:
            cluster_names_to_genomes[bgc] = genomes.index(identifier)
    run.run_data["input"]["accession"] = [{"id": "genome_{}".format(i), "label": acc} for i, acc in enumerate(genomes)]
    run.run_data["input"]["accession_newick"] = [] # todo ...
    run.run_data["input"]["classes"] = [{"label": cl} for cl in classes] # todo : colors
    run.run_data["input"]["bgc"] = [{"id": cluster_names[idx], "acc": cluster_names_to_genomes[cluster_names[idx]], "class": cluster_names_to_classes[cluster_names[idx]]} for idx in input_clusters_idx]

def update_family_data(RUNDATA_NETWORKS_PER_RUN, INPUT_CLUSTERS_IDX, CLUSTER_NAMES, MIBIG_SET):
    for network_key in RUNDATA_NETWORKS_PER_RUN:
        for network in RUNDATA_NETWORKS_PER_RUN[network_key]:
            for family in network["families"]:
                new_members = []
                mibig = []
                for bgcIdx in family["members"]:
                    if bgcIdx in INPUT_CLUSTERS_IDX:
                        new_members.append(INPUT_CLUSTERS_IDX.index(bgcIdx))
                    else: # is a mibig bgc
                        clusterName = CLUSTER_NAMES[bgcIdx]
                        if clusterName in MIBIG_SET:
                            mibig.append(clusterName)
                family["mibig"] = mibig
                family["members"] = new_members


def generate_results_per_cutoff_value(run, rundata_networks_per_run, html_subs_per_run):
    for cutoff in run.cluster.cutoff_list:
        # update overview.html
        html_folder_for_this_cutoff = "{}_c{:.2f}".format(run.directories.network_html, cutoff)
        run_data_for_this_cutoff = run.run_data.copy()
        run_data_for_this_cutoff["networks"] = rundata_networks_per_run[html_folder_for_this_cutoff]
        with open(os.path.join(html_folder_for_this_cutoff, "run_data.js"), "w") as run_data_js:
            run_data_js.write("var run_data={};\n".format(json.dumps(run_data_for_this_cutoff, indent=4, separators=(',', ':'), sort_keys=True)))
            run_data_js.write("dataLoaded();\n")
        # update bgc_results.js
        RUN_STRING = "{}_c{:.2f}".format(run.run_name, cutoff)
        RESULTS_PATH = os.path.join(run.directories.output, "html_content", "js",
                                    "bigscape_results.js")
        add_to_bigscape_results_js(RUN_STRING, html_subs_per_run[html_folder_for_this_cutoff],
                                      RESULTS_PATH)
