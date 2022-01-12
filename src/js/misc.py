import os
import json


def add_to_bigscape_results_js(module_name, subs, result_js_file):
    bigscape_results = []
    if os.path.isfile(result_js_file):
        with open(result_js_file, "r") as bs_js:
            line = bs_js.read()
            assert line.startswith("var bigscape_results = ")
            assert line.endswith(";")
            bigscape_results = json.loads(line[23:-1])
    bigscape_results.append({"label" : module_name, "networks" : subs})
    with open(result_js_file, "w") as bs_js:
        json_string = json.dumps(bigscape_results, indent=4, separators=(',', ':'), sort_keys=True)
        bs_js.write("var bigscape_results = {};".format(json_string))
