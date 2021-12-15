import os
from collections import defaultdict

def stockholm_parser(stkFile):
    reference = ""
    algnDict = defaultdict(str)
    algnFile = stkFile[:-3] + 'algn'
    if not os.path.isfile(algnFile): # prevents overwriting algn files
        with open(stkFile, 'r') as infile:
            for l in infile:
                line = l.strip()
                if line.startswith("#=GC RF"):
                    reference += line[7:].strip()
                elif line == "":
                    continue
                elif line[0] == "/" or line[0] == "#":
                    continue
                else:
                    a = line.split(" ")
                    header = a[0]
                    algn = a[-1]
                    algnDict[header] += algn
                    
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
                slicing_tuples.append((start,pos))
        if state_reference:
            slicing_tuples.append((start, len(reference)))
    
    if len(algnDict) > 0:
        with open(algnFile, "w") as outfile:
            for header in algnDict:
                sequence = ""
                for a,b in slicing_tuples:
                    sequence += algnDict[header][a:b]
                outfile.write(">{}\n".format(header))
                outfile.write(sequence + "\n")
    return
