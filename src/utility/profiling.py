"""Module to contain profiling code.

Authors:
    Arjan Draisma (arjan.draisma@wur.nl)
"""

from multiprocessing import Process, Queue
import os
import threading
from time import sleep
import time
import psutil

def get_stats(process: psutil.Process):
    """Method to return a set of current usage statistics for a given process
    """
    cpu = process.cpu_percent()
    mem_mb = process.memory_info().vms / 1000000
    mem_percent = process.memory_percent()
    return cpu, mem_mb, mem_percent


def collect_consumption(log_path: str, command_queue: Queue, update_interval: int):
    """Worker thread to periodically report cpu and memory usage"""
    with open(log_path, "w", encoding="UTF-8") as profile_log:
        profile_log.write("time,cpu,processes,mem_used_mb,mem_used_perc\n")
        while threading.main_thread().is_alive():
            prefix = time.strftime("%Y-%m-%d %H:%M:%S:%sss", time.localtime())
            if not command_queue.empty():
                command, args = command_queue.get()
                if command == 1:
                    break

            # main thread
            main_process = psutil.Process(os.getpid()).parent()
            cpu, mem_mb, mem_percent = get_stats(main_process)

            # children
            for c_process in main_process.children(recursive=True):
                c_cpu, c_mem_mb, c_mem_perc = get_stats(c_process)
                mem_mb += c_mem_mb
                mem_percent += c_mem_perc
            
            processes = len(main_process.children(recursive=True))

            cpu = main_process.cpu_percent(update_interval)

            # format and write log line
            log_line = f"{prefix},{cpu},{processes},{mem_mb:.2f},{mem_percent:.4f}\n"
            profile_log.write(log_line)
            profile_log.flush()

class Profiler:
    """Class to wrap the profiler functionality"""
    worker: Process
    command_queue: Queue


    def __init__(self, options, log_timestamp):
        
        # look for override
        if options.log_path is None:
            log_path = options.outputdir
        else:
            log_path = options.log_path

        log_file = os.path.join(log_path, log_timestamp + "_profile.log")
        self.command_queue = Queue()
        self.worker = Process(target=collect_consumption, args=(log_file, self.command_queue, 0.5))

    def start(self):
        """Starts the worker thread"""
        self.worker.start()

    def stop(self):
        """Stops the worker thread"""
        self.command_queue.put((1, None))
        self.worker.join()
