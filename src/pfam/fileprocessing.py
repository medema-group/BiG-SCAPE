import os
import json

from src.pfam.misc import generatePfamColorsMatrix

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

def create_pfam_js(run, pfam_info):
    PFAMS_JS_FILE = os.path.join(run.directories.output, "html_content", "js", "pfams.js")
    if not os.path.isfile(PFAMS_JS_FILE):
        with open(PFAMS_JS_FILE, "w") as pfams_js:
            PFAM_JSON = {}
            PFAM_COLORS = generatePfamColorsMatrix(os.path.join(os.path.dirname(os.path.realpath(__file__)), "domains_color_file.tsv"))
            for pfam_code in pfam_info:
                pfam_obj = {}
                if pfam_code in PFAM_COLORS:
                    pfam_obj["col"] = PFAM_COLORS[pfam_code]
                else:
                    pfam_obj["col"] = "255,255,255"
                pfam_obj["desc"] = pfam_info[pfam_code][1]
                PFAM_JSON[pfam_code] = pfam_obj
            pfams_js.write("var pfams={};\n".format(json.dumps(PFAM_JSON, indent=4, separators=(',', ':'), sort_keys=True)))