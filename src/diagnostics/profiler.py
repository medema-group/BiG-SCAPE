"""Module containing code to handle profiling"""

# from python
import threading
import os
import time
from datetime import datetime
from pathlib import Path
from multiprocessing import Process, Queue
import psutil


class Profiler:
    """
    Class containing the profiler functionality

    Attributes:
        worker: Process
        command_queue: Queue
    """

    def __init__(self, profile_path: Path) -> None:
        self.command_queue: Queue = Queue()
        self.worker = Process(
            target=collect_consumption, args=(profile_path, self.command_queue, 0.5)
        )

    def start(self) -> None:
        """Starts the worker thread"""
        self.worker.start()

    def stop(self):
        """Stops the worker thread"""
        self.command_queue.put((1, None))
        self.worker.join()


def calc_cpu_percent(
    process: psutil.Process, start_time: float, update_interval: float
):
    """Calculate cpu percent by measuring cpu_time over elapsed time,
        i.e. difference between the cpu time at the start and end of the interval,
        and convert to percentage by dividing over elapsed interval

    Args:
        process (psutil.Process): any process
        start_time (float): CPU time at the start of interval, in seconds
        end_time (float): CPU time at the end of interval, in seconds
        update_interval (float): in seconds

    Returns:
        float: cpu percent usage
    """

    end_time = process.cpu_times().user
    elapsed_time = end_time - start_time
    cpu_percent = elapsed_time / update_interval

    return cpu_percent


def get_stats(process: psutil.Process, start_time: float, update_interval: float):
    """Return a set of usage stats for a given process

    Args:
        process (psutil.Process): any process
        start_time (float): CPU time at the start of interval, in seconds
        end_time (float): CPU time at the end of interval, in seconds
        update_interval (float): in seconds

    Returns:
        tuple: cpu_percent, mem_mb, mem_percent
    """

    mem_mb = process.memory_info().rss / 1000000
    mem_percent = process.memory_percent()

    cpu_percent = calc_cpu_percent(process, start_time, update_interval)

    return (cpu_percent, mem_mb, mem_percent)


def child_start_cpu_time(main_process: psutil.Process):
    """Generate dictionary CPU times of the child processes
        at the start of each interval

    Args:
        main_process (psutil.Process): main process

    Returns:
        dict: {child_id: user cpu_time in seconds}
    """

    set_child_cpu_time_dict = {}
    for child_process in main_process.children(recursive=True):
        c_pid = f"CHILD_{child_process.pid}"
        set_child_cpu_time_dict[c_pid] = child_process.cpu_times().user

    return set_child_cpu_time_dict


def collect_consumption(
    profile_path: Path, command_queue: Queue, update_interval: int
) -> None:
    """Worker thread to periodically report cpu and memory usage"""
    mem_dict: dict = {}
    cpu_dict: dict = {}

    with open(profile_path, "w", encoding="UTF-8") as profile_log:
        profile_log.write("Time,Type,CPU,Processes,mem_used_mb,mem_used_perc\n")

        start_time = time.time()

        main_process = psutil.Process(os.getpid()).parent()

        while threading.main_thread().is_alive():
            if not command_queue.empty():
                command, args = command_queue.get()
                if command == 1:
                    break

            # cpu time at the start of the interval for main process
            set_main_cpu_time = main_process.cpu_times().user
            # and for children processes
            set_child_cpu_time_dict = child_start_cpu_time(main_process)

            # elapse time interval
            time.sleep(update_interval)

            # now that the interval has passed, we can get the usage stats:
            # cpu_percent is measured via the cpu_time of this interval
            # and memory is measured at the timepoint where the interval ends

            # set timestamp to be used in human readable report
            prefix_time = datetime.now().isoformat(sep=" ", timespec="milliseconds")
            # set time stamp to be used in graph reports
            current_time = time.time()
            elapsed_time = current_time - start_time

            # start dicts to record stats per time stamp to be used in graph reports
            mem_dict[elapsed_time] = {}
            cpu_dict[elapsed_time] = {}

            # get main process usage stats
            m_cpu_percent, m_mem_mb, m_mem_percent = get_stats(
                main_process, set_main_cpu_time, update_interval
            )

            # populate stats to be used in graph reports
            mem_dict[elapsed_time]["MAIN"] = m_mem_mb
            cpu_dict[elapsed_time]["MAIN"] = m_cpu_percent

            # print 'MAIN' line to profile report
            main_line = f"{prefix_time},MAIN,{m_cpu_percent},1,{m_mem_mb:.2f},{m_mem_percent:.2f}\n"
            profile_log.write(main_line)

            # start multi/cumulative stats
            total_cpu_percent = m_cpu_percent
            total_mem_mb = m_mem_mb
            total_mem_percent = m_mem_percent
            child_process_count = len(main_process.children(recursive=True))

            for child_process in main_process.children(recursive=True):
                c_pid = f"CHILD_{child_process.pid}"
                try:
                    # in case a child process was created between the time of generating the set_child_cpu_time_dict
                    # and now, i.e. during the elapsed interval
                    if c_pid not in set_child_cpu_time_dict:
                        continue

                    # get child process usage stats
                    c_cpu_percent, c_mem_mb, c_mem_percent = get_stats(
                        child_process, set_child_cpu_time_dict[c_pid], update_interval
                    )

                    # populate stats to be used in graph reports
                    mem_dict[elapsed_time][c_pid] = c_mem_mb
                    cpu_dict[elapsed_time][c_pid] = c_cpu_percent

                    # write CHILD line to profile report
                    child_line = (
                        f"{prefix_time},{c_pid},{c_cpu_percent},1,{c_mem_mb:.2f},"
                        f"{c_mem_percent:.2f}\n"
                    )
                    profile_log.write(child_line)

                except psutil.NoSuchProcess:
                    continue

                # update multi/cumulative stats
                total_cpu_percent += c_cpu_percent
                total_mem_mb += c_mem_mb
                total_mem_percent += c_mem_percent

            # write MULTI line to profile report
            multi_line = (
                f"{prefix_time},MULTI,{total_cpu_percent},{child_process_count},"
                f"{total_mem_mb:.2f},{m_mem_percent:.2f}\n"
            )
            profile_log.write(multi_line)
            profile_log.flush()
