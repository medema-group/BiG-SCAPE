import os
from collections import defaultdict

def stockholm_parser(stk_file):
    """shoo"""
    reference = ""
    algn_dict = defaultdict(str)
    algn_file = stk_file[:-3] + 'algn'
    if not os.path.isfile(algn_file): # prevents overwriting algn files
        with open(stk_file, 'r') as infile:
            for line in infile:
                line = line.strip()
                if line.startswith("#=GC RF"):
                    reference += line[7:].strip()
                elif line == "":
                    continue
                elif line[0] == "/" or line[0] == "#":
                    continue
                else:
                    start = line.split(" ")
                    header = start[0]
                    algn = start[-1]
                    algn_dict[header] += algn

        # get start-end coordinates of every "x" island (original consensus)
        # in the reference
        state_reference = False
        slicing_tuples = []
        for pos in range(len(reference)):
            if reference[pos] == "x" and not state_reference:
                state_reference = True
                start = pos
            if reference[pos] == "." and state_reference:
                state_reference = False
                end = pos
                slicing_tuples.append((start, end))
        if state_reference:
            slicing_tuples.append((start, len(reference)))

    if len(algn_dict) > 0:
        with open(algn_file, "w") as outfile:
            for header in algn_dict:
                sequence = ""
                for start, end in slicing_tuples:
                    sequence += algn_dict[header][start:end]
                outfile.write(">{}\n".format(header))
                outfile.write(sequence + "\n")
    return
