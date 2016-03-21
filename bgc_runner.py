#!/usr/bin/env python


"""
Student/programmer: Marley Yeong
marleyyeong@live.nl
supervisor: Marnix Medema

Usage: run bgc_networks.py with different parameter combinations

"""
from optparse import OptionParser
from subprocess import Popen
from itertools import islice
import os
from functions import frange

def param_combinations(outputdir, steps):
    commands = []
    for i in frange(0,1,steps):

        for j in frange(0,1,steps):

            for k in frange(0,1,steps):

                for l in frange(0,1,steps):
                    if i+j+k == 1:
                        cmd = "python ~/bgc_networks/bgc_networks.py -o " + outputdir + " --sim_cutoffs \"0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75\" \
                        --skip_hmmscan --Jaccardw " + str(i) + " --DDSw " + str(j) +" --GKw " + str(k) + " -a " + str(l) + " -s --al_method \"--retree 2\" --maxiterate 10 --include_disc_nodes"
                        commands.append(cmd)
    
                
    max_workers = 8  # no more than 8 concurrent processes
    processes = (Popen(cmd, shell=True) for cmd in commands)
    running_processes = list(islice(processes, max_workers))  # start new processes
    while running_processes:
        for i, process in enumerate(running_processes):
            if process.poll() is not None:  # the process has finished
                try:
                    running_processes[i] = next(processes)  # start new process
                except StopIteration:
                    del running_processes[i]
                    break
                

#===============================================================================
# def param_combinations(outputdir, steps):
# 
#     for i in bgc_functions.frange(0,1,steps):
# 
#         for j in bgc_functions.frange(0,1,steps):
# 
#             for k in bgc_functions.frange(0,1,steps):
# 
#                 for l in bgc_functions.frange(0,1,steps):
#                     if i+j+k == 1:
#                         cmd = "python ~/bgc_networks/bgc_networks.py -o " + outputdir + " --sim_cutoffs \"0.4\" \
#                         --skip_hmmscan --Jaccardw " + str(i) + " --DDSw " + str(j) +" --GKw " + str(k) + " -a " + str(l) + " -s"
#     
#                 
#                         #print cmd
#                         processes = set() #Will only contain one subprocess at a time. Probably because
#                                           #the docker container hides the process from subprocess, so it seems
#                                           #as if the process is immediately finished. This variable is still needed,
#                                           #otherwise the function will not wait for the last antismash instance to be finished.
#                                           #It is possible that more than <max_processes> are executed if the first files from 
#                                           #the batch are much bigger than the last file in the batch.
#                                             
#                         processes.add(subprocess.Popen(cmd, shell=True))
#             
#     #Check if all the child processes were closed
#     for p in processes:
#         if p.poll() is None:
#             p.wait()  
#                 
#===============================================================================


def CMD_parser():
    parser = OptionParser()
    parser.add_option("-n", "--networkfile", dest="networkfile", default="",
                      help="name of networkfile")
    parser.add_option("-o", "--outputdir", dest="outputdir", default="second_test",
                      help="output directory, this contains your pfd,pfs,network and hmmscan output files")
    parser.add_option("--steps", dest="steps", default=0.1,
                      help="DDS weight")
 
    (options, args) = parser.parse_args()
    return options, args
 
 
 
if __name__=="__main__":
     
    options, args = CMD_parser()
    param_combinations(options.outputdir, float(options.steps))
    #===========================================================================
    # networkfile = options.networkfile
    # groups = network_reader.load_groups("groups.txt")
    # print network_reader.get_network_score(groups, networkfile)
    #===========================================================================
    #param_combinations(0.2)
    
