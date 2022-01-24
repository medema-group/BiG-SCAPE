from multiprocessing.dummy import Array
import os
import sys

import src.bgctools as bgctools
import src.utility as utility

# disable pyling complaining about backwards compatibility measures
#pylint: disable=wrong-import-order,redefined-builtin,invalid-name,undefined-variable,import-error
from sys import version_info
if version_info[0] == 2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0] == 3:
    import pickle # for storing and retrieving dictionaries
#pylint: enable=wrong-import-order,redefined-builtin,invalid-name,undefined-variable,import-error


class BgcInfo:
    # dictionary of this structure:
    # BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
    #  'specific_domain_name_2'] } }
    # - cluster_name_x: cluster name (can be anything)
    # - general_domain_name_x: PFAM ID, for example 'PF00550'
    # - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names
    # in DMS unequivocally. e.g. 'PF00550_start_end', where start and end are genomic positions
    domain_name_info = {}

    def load_from_file(self, run):
        try:
            with open(os.path.join(run.directories.cache, "BGCs.dict"), "r") as bgc_file:
                self.domain_name_info = pickle.load(bgc_file)
                bgc_file.close()
        except IOError:
            sys.exit("BGCs file not found...")

    def load_pfds(self, run, cluster_base_names, try_resume_multiple_alignment):
        for outputbase in cluster_base_names:
            if run.options.verbose:
                print("   Processing: " + outputbase)

            pfd_file = os.path.join(run.directories.pfd, outputbase + ".pfd")
            filtered_matrix = [[part.strip() for part in l.split('\t')] for l in open(pfd_file)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(run.directories.bgc_fasta, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not try_resume_multiple_alignment:
                with open(fasta_file, "r") as fasta_file_handle:
                    # all fasta info from a BGC
                    fasta_dict = utility.fasta_parser(fasta_file_handle)
                utility.save_domain_seqs(filtered_matrix, fasta_dict,
                                         run.directories.domains, outputbase)

            self.domain_name_info[outputbase] = bgctools.generate_domain_name_info_dict(filtered_matrix)

            del filtered_matrix[:]

    def save_to_file(self, run):
        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(run.directories.cache, "BGCs.dict"), "wb") as bgc_file:
            pickle.dump(self.domain_name_info, bgc_file)
            bgc_file.close()
