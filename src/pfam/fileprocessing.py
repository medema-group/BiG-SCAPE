import os

def parse_pfam_a(run):
    pfam_info = {}
    with open(os.path.join(run.directories.pfam, "Pfam-A.hmm"), "r") as pfam:
        put_in_dict = False
        # assuming that the order of the information never changes
        for line in pfam:
            if line[:4] == "NAME":
                name = line.strip()[6:]
            if line[:3] == "ACC":
                acc = line.strip()[6:].split(".")[0]
            if line[:4] == "DESC":
                desc = line.strip()[6:]
                put_in_dict = True

            if put_in_dict:
                put_in_dict = False
                pfam_info[acc] = (name, desc)
    return pfam_info
