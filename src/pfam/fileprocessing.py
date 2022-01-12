import os
import json

from src.pfam.misc import generate_pfam_colors_matrix

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
    pfams_js_file = os.path.join(run.directories.output, "html_content", "js", "pfams.js")
    if not os.path.isfile(pfams_js_file):
        with open(pfams_js_file, "w") as pfams_js:
            pfam_json = {}
            domains_colors_path = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                               "domains_color_file.tsv")
            pfam_colors = generate_pfam_colors_matrix(domains_colors_path)
            for pfam_code in pfam_info:
                pfam_obj = {}
                if pfam_code in pfam_colors:
                    pfam_obj["col"] = pfam_colors[pfam_code]
                else:
                    pfam_obj["col"] = "255,255,255"
                pfam_obj["desc"] = pfam_info[pfam_code][1]
                pfam_json[pfam_code] = pfam_obj
            json_string = json.dumps(pfam_json, indent=4, separators=(',', ':'), sort_keys=True)
            pfams_js.write("var pfams={};\n".format(json_string))
