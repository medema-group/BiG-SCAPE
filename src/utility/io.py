"""Module containing various input/output utilities

Authors: Jorge Navarro, Arjan Draisma
"""


import os
import sys
import logging

def create_directory(path, kind, clean):
    """Create a directory at a given path
    
    parameters:
        path: path to create
        kind: description of the kind of folder that is being created
        clean: whether or not to remove the contents of the folder if it
        already exists
    """
    # TODO consider makedirs(path,exists_ok=True)
    try:
        os.makedirs(path)
    except OSError as exception:
        # 17 (Linux): "[Errno 17] File exists";
        # 183 (Windows) "[Error 183] Cannot create a file when that file already exists"
        if "Errno 17" in str(exception) or "Error 183" in str(exception):
            logging.info(" %s folder already exists", kind)
            if clean:
                logging.info("  Cleaning folder")
                for thing in os.listdir(path):
                    os.remove(os.path.join(path, thing))
        else:
            logging.error("unexpected error creating " + kind + " folder")
            logging.error(exception.strerror)
            sys.exit(1)

def write_parameters(output_folder, parameters):
    """Write a file with all the details of the run.
    Will overwrite previous versions"""

    with open(os.path.join(output_folder, "parameters.txt"), "w") as parameters_file:
        parameters_file.write(" ".join(parameters))
