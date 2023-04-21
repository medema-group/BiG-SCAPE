"""Module containing code to handle profiling"""

# from python
from pathlib import Path
from multiprocessing import Process, Queue
import psutil


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


def collect_consumption() -> None:
    return None


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
