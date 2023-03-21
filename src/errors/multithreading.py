"""Contains errors for the multithreading part of BiG-SCAPE"""


class WorkerPoolSetupError(Exception):
    """Thrown when no function is assigned to a worker pool"""

    def __init__(self, *args: object) -> None:
        super().__init__(*args)


class WorkerSetupError(Exception):
    """Thrown when operations are performed on an unstarted worker"""

    def __init__(self, *args: object) -> None:
        super().__init__(*args)
