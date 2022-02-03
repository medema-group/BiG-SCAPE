import os
import time
import logging

def timeit(funct):
    """Writes the runtimes of functions to a file and prints them on screen.

    Input:
    - funct: a function
    Output:
    - That function its output.
    - Runtime of that function in a file called commands.txt and on screen.
    """
    def _wrap(*args):
        start_time = time.time()
        ret = funct(*args)
        runtime = time.time()-start_time
        runtime_string = '{} took {:.3f} seconds'.format(funct.__name__, runtime)
        with open(os.path.join(log_folder, "runtimes.txt"), 'a') as timings_file:
            timings_file.write(runtime_string + "\n")
        logging.info(runtime_string)
        return ret
    
    return _wrap