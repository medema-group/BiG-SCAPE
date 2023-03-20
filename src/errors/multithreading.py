"""Contains errors for the multithreading part of BiG-SCAPE"""


class NoFunctionError(Exception):
    """Thrown when no function is assigned to a worker pool"""

    def __init__(self) -> None:
        super().__init__("No function assigned to worker pool!")


class WorkerNotStartedError(Exception):
    """Thrown when operations are performed on an unstarted worker"""

    def __init__(self, id) -> None:
        super().__init__("Worker with id %d process was not started!", id)
