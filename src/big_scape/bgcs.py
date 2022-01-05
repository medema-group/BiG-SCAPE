import os
import sys

import src.bgctools as bgctools
import src.utility as utility

from sys import version_info
if version_info[0] == 2:
    range = xrange
    import cPickle as pickle # for storing and retrieving dictionaries
elif version_info[0] == 3:
    import pickle # for storing and retrieving dictionaries

class BGCS:
    # dictionary of this structure:
    # BGCs = {'cluster_name_x': { 'general_domain_name_x' : ['specific_domain_name_1',
    #  'specific_domain_name_2'] } }
    # - cluster_name_x: cluster name (can be anything)
    # - general_domain_name_x: PFAM ID, for example 'PF00550'
    # - specific_domain_name_x: ID of a specific domain that will allow to you to map it to names
    # in DMS unequivocally. e.g. 'PF00550_start_end', where start and end are genomic positions
    bgc_dict = {}

    def load_from_file(self, RUN):
        try:
            with open(os.path.join(RUN.directories.cache, "BGCs.dict"), "r") as BGC_file:
                self.bgc_dict = pickle.load(BGC_file)
                BGC_file.close()
        except IOError:
            sys.exit("BGCs file not found...")
        
    def load_pfds(self, RUN, CLUSTER_BASE_NAMES, try_resume_multiple_alignment):
        for outputbase in CLUSTER_BASE_NAMES:
            if RUN.options.verbose:
                print("   Processing: " + outputbase)

            pfdFile = os.path.join(RUN.directories.pfd, outputbase + ".pfd")
            FILTERED_MATRIX = [[part.strip() for part in l.split('\t')] for l in open(pfdFile)]

            # save each domain sequence from a single BGC in its corresponding file
            fasta_file = os.path.join(RUN.directories.bgc_fasta, outputbase + ".fasta")

            # only create domain fasta if the pfd content is different from original and
            #  domains folder has been emptied. Else, if trying to resume alignment phase,
            #  domain fasta files will contain duplicate sequence labels
            if not try_resume_multiple_alignment:
                with open(fasta_file, "r") as fasta_file_handle:
                    # all fasta info from a BGC
                    fasta_dict = utility.fasta_parser(fasta_file_handle)
                utility.save_domain_seqs(FILTERED_MATRIX, fasta_dict,
                                         RUN.directories.domains, outputbase)

            self.bgc_dict[outputbase] = bgctools.bgc_dict_gen(FILTERED_MATRIX)

            del FILTERED_MATRIX[:]

    def save_to_file(self, RUN):
        # store processed BGCs dictionary for future re-runs
        with open(os.path.join(RUN.directories.cache, "BGCs.dict"), "wb") as BGC_file:
            pickle.dump(self.bgc_dict, BGC_file)
            BGC_file.close()
