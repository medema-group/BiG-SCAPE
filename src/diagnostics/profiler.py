"""Module containing code to handle profiling"""

# from python
from pathlib import Path
from multiprocessing import Process, Queue


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
