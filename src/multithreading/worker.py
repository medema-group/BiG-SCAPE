"""Class to contain methods to create and destroy pools of workers for multithreading
tasks
"""

# from python
from __future__ import annotations
import logging
from multiprocessing import Pipe, Process, cpu_count
from multiprocessing.connection import Connection, wait
from typing import Any, Optional, Callable, cast
from random import choice

# from other modules
from src.errors import WorkerPoolSetupError, WorkerSetupError


class WorkerPool:
    """Helper class to contain a pool of workers, used handle large workloads

    The most important aspect of this class is the self.worker_function, which is the
    function that is passed to worker objects. This function only ever has one argument
    (a connection object created by the multiprocessing.pipe command).

    This connection should only expect tuples to come in. the first element of a tuple
    will be an int that specifies the command type:

    0. STOP - stop this tread please. no other elements in tuple
    1. DEGUB - show a debug message. the second element will be the message
    2. TASK - perform a task. what the rest of this tuple looks like is up to the
    individual implementation

    Connections can also receive data back from the worker. Possible responses are:

    3. RESULT - Receive a result from a task. Again the rest of the elements in this
    tuple are up to individual implementation
    4. STOP - This process is stopping
    5. ERROR - An error has occured in this process. second element should have some
    message
    """

    def __init__(self, num_workers=cpu_count()) -> None:
        self.workers = [Worker(worker_id, self) for worker_id in range(num_workers)]

        # list of connections, which we will monitor for traffic from the workers
        self.connections: list[Connection] = []

        # the function we will pass to the worker when it is created
        self.worker_function: Optional[Callable] = None

    def start(self) -> None:
        if self.worker_function is None:
            raise WorkerPoolSetupError(
                "Tried to start worker pool without a worker function"
            )

        for worker in self.workers:
            logging.debug("Starting worker with id %d", worker.id)
            pool_connection, worker_connection = Pipe(duplex=True)
            self.connections.append(pool_connection)
            worker.start(pool_connection, worker_connection)

            task_data = (worker.id,)
            worker.send_task(Worker.START, task_data)

    def stop(self) -> None:
        """Stop the worker pool, sending a stop signal to all workers"""
        for worker in self.workers:
            worker.stop()

    def kill(self) -> None:
        """Force-stops the workers. used in clean up if an error occurs"""
        for worker in self.workers:
            worker.stop(True)

    def join(self, timeout: Optional[float] = None) -> None:
        """calls .join() on all child workers, waiting for the processes to close"""
        for worker in self.workers:
            worker.join(timeout)

    def send_task_one(self, data: tuple[Any]):
        """Sends a task to a random worker"""
        connection = choice(self.connections)
        request = (Worker.TASK,) + data
        connection.send(request)

    def read_one(self, timeout: Optional[float] = None):
        available_connections = cast(list[Connection], wait(self.connections, timeout))
        connection = choice(available_connections)
        return connection.recv()

    def wait(self) -> list[Connection]:
        return cast(list[Connection], wait(self.connections))


class Worker:
    """Class that takes a function and takes data from a pool to execute on the function"""

    STOP = 0  # command to stop the process, or a notification that a worker stopped
    START = 1  # command to ask the function to print a debug
    TASK = 2  # command to perform a task
    RESULT = 3  # response with a result
    ERROR = 5  # response announcing an error occurred

    def __init__(self, id: int, pool: WorkerPool):
        self.id = id
        self.pool = pool

        self.process: Optional[Process] = None

        self.connection: Optional[Connection] = None

    def start(self, pool_connection: Connection, worker_connection: Connection):
        """Starts the worker, creating a process and starting the listening for input

        Args:
            connection (Connection): Connection for the worker created by the pool
        """
        logging.debug("Starting process for worker with id %d", self.id)

        self.process = Process(
            target=self.pool.worker_function,
            args=(
                self.id,
                worker_connection,
            ),
        )
        self.connection = pool_connection

        self.process.start()

    def stop(self, force=False):
        """Stops the process in this worker. if force=True, kills the process without
        sending a nice stop message
        """
        if force:
            logging.warning("Force-closing process for worker %d", self.id)
            self.process.kill()

        logging.debug("Sending stop signal to worker %d", self.id)
        self.connection.send((Worker.STOP, self.id))

    def join(self, timeout: Optional[float] = None):
        """Wrapper to call .join on this worker process

        Args:
            timeout (float, optional): Time in seconds to wait for process to close.
            Defaults to None.
        """
        if self.process is None:
            raise WorkerSetupError(
                "Worker with id %d process was not started!", self.id
            )

        self.process.join(timeout)

    def send_task(self, command: int, data: tuple[Any]):
        """Sends a command to the worker process

        Args:
            command (int): The command to send. E.g. Worker.TASK
        """

        if self.connection is None:
            raise WorkerSetupError(
                "Worker with id %d process was not started!", self.id
            )

        logging.debug("Sending message to worker %d", self.id)

        # convert input arguments to a tuple
        command_package = (command,) + data

        self.connection.send(command_package)
