"""Class to contain methods to create and destroy pools of workers for multithreading
tasks
"""

# from python
import logging
from multiprocessing import Pipe, Process, cpu_count, Queue
from multiprocessing.connection import Connection
from typing import Optional, Callable

# from other modules
from src.errors import NoFunctionError, WorkerNotStartedError


class WorkerPool:
    """Helper class to contain a pool of workers, used to queue tasks and handle large
    workloads

    The most important aspect of this class is the self.worker_function, which is the
    function that is passed to worker objects. This function only ever has one argument
    (a connection object created by the multiprocessing.pipe command).

    This connection should only expect tuples to come in. Furthermore, this function
    should at least expect these two commands to come in at any point through the
    command:

    ("_STOP") - stops the current thread

    ("_ANNOUNCE", id) - logs "Process for worker {id} has started" to debug

    ("_TASK", **argv) - an actual task for the worker function to execute
    """

    def __init__(self, num_workers=cpu_count()):
        self.workers = [Worker(worker_id, self) for worker_id in range(num_workers)]

        # queues for tasks
        self.input_queue = Queue()
        self.output_queue = Queue()

        # list of connections, which we will monitor for traffic from the workers
        self.connections: list[Connection] = []

        # the function we will pass to the worker when it is created
        self.worker_function: Optional[Callable] = None

    def queue_task(self, arguments: tuple):
        """Adds a task to the queue for this pool of workers. Queues work on a FIFO
        basis

        Args:
            arguments (tuple): A tuple of arguments to pass to worker_function. The
            formatting of the tuple needs to be able to be processed by the worker
        """
        self.input_queue.put(arguments)

    def start(self) -> None:
        if self.worker_function is None:
            raise NoFunctionError()

        for worker in self.workers:
            logging.debug("Starting worker with id %d", worker.id)
            worker.start()

    def stop(self) -> None:
        """Stop the worker pool, sending a kill signal to all workers"""
        for worker in self.workers:
            worker.stop()

    def join(self, timeout: Optional[float] = None) -> None:
        """calls .join() on all child workers, waiting for the processes to close"""
        for worker in self.workers:
            worker.join(timeout)


class Worker:
    """Class that takes a function and takes data from a pool to execute on the function"""

    COMM_STOP = "_STOP"  # command to stop the process
    COMM_DEBUG = "_DEBUG"  # command to ask the function to print a debug
    COMM_DEBUG = "_TASK"  # command to perform a task

    def __init__(self, id: int, pool: WorkerPool):
        self.id = id
        self.pool = pool

        self.process: Optional[Process] = None

        self.connection: Optional[Connection] = None

    def start(self):
        """Starts the worker, creating a process and starting the listening for input"""
        # create a new pipe. this consists of two connections: one for the pool and one
        # for the worker. either side can send and receive stuff through its connection
        logging.debug("Starting process for worker with id %d", self.id)
        pool_connection, worker_connection = Pipe(duplex=True)

        self.process = Process(
            target=self.pool.worker_function, args=(worker_connection,)
        )

        self.process.start()
        self.connection = pool_connection

        # as described in the worker_pool class, we expect the function worker to be

    def stop(self, force=False):
        """Stops the process in this worker. if force=True, kills the process without
        sending a nice kill message
        """
        if force:
            logging.debug("Force-closing process for worker %d", self.id)
            self.process.kill()

        logging.debug("Sending stop signal to worker %d", self.id)
        self.send(Worker.COMM_STOP)

    def join(self, timeout: Optional[float] = None):
        """Wrapper to call .join on this worker process

        Args:
            timeout (float, optional): Time in seconds to wait for prceoss to close.
            Defaults to None.
        """
        if self.process is None:
            raise WorkerNotStartedError(self.id)

        self.process.join(timeout)

    def send(self, command: str, **argv):
        """Sends a command to the worker process

        Args:
            command (str): The command to send. E.g. Worker.COMM_STOP
        """

        if self.connection is None:
            raise WorkerNotStartedError(self.id)

        # convert input arguments to a tuple
        command_package = (command,) + tuple(argv.items())

        self.connection.send(command_package)
