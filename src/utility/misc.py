import logging


def get_anchor_domains(filename):
    """Get the anchor/marker domains from a txt file.
    This text file should contain one Pfam id per line.
    A second column (separated by a tab) with more comments is allowed"""

    domains = set()

    try:
        with open(filename, "r") as handle:
            for line in handle:
                # handle comments and empty lines
                if line[0] != "#" and line.strip():
                    # ignore domain versions
                    domains.add(line.strip().split("\t")[0].split(".")[0])
        return domains
    except IOError:
        logging.warning("You have not provided the anchor_domains.txt file.")
        logging.warning("if you want to make use of the anchor domains in the DSS distance\
            metric, make a file that contains a Pfam domain on each line.")
        return set()
