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


def collect_consumption(log_path: str, command_queue: Queue, update_interval: int):
    """Worker thread to periodically report cpu and memory usage"""
    with open(log_path, "w", encoding="UTF-8") as profile_log:
        profile_log.write("# Note: this displays total usage on the system, and may be affected \
                          by other processes running on the machine.")
        profile_log.write("time,cpu,mem_used_mb,mem_used_perc\n")
        while threading.main_thread().is_alive():
            prefix = time.strftime("%Y-%m-%d_%H-%M-%S", time.localtime())
            if not command_queue.empty():
                command, args = command_queue.get()
                if command == 1:
                    break

            cpu = psutil.cpu_percent()
            mem = psutil.virtual_memory()
            profile_log.write(f"{prefix},{cpu},{mem.used / 1000000:.2f},{mem.percent}\n")

            sleep(update_interval)

class Profiler:
    """Class to wrap the profiler functionality"""
    worker: Process
    command_queue: Queue


    def __init__(self, options):
        log_file = os.path.join(options.log_path, "profile.log")
        self.command_queue = Queue()
        self.worker = Process(target=collect_consumption, args=(log_file, self.command_queue, 1))

    def start(self):
        """Starts the worker thread"""
        self.worker.start()

    def stop(self):
        """Stops the worker thread"""
        self.command_queue.put((1, None))
        self.worker.join()
