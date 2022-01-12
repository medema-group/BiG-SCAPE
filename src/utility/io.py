import os
import sys

def create_directory(path, kind, clean):
    # TODO consider makedirs(path,exists_ok=True)
    try:
        os.makedirs(path)
    except OSError as exception:
        # 17 (Linux): "[Errno 17] File exists";
        # 183 (Windows) "[Error 183] Cannot create a file when that file already exists"
        if "Errno 17" in str(exception) or "Error 183" in str(exception):
            print(" " + kind + " folder already exists")
            if clean:
                print("  Cleaning folder")
                for thing in os.listdir(path):
                    os.remove(os.path.join(path, thing))
        else:
            print("Error: unexpected error creating " + kind + " folder")
            sys.exit(str(exception))

def write_parameters(output_folder, parameters):
    """Write a file with all the details of the run.
    Will overwrite previous versions"""

    with open(os.path.join(output_folder, "parameters.txt"), "w") as parameters_file:
        parameters_file.write(" ".join(parameters))
