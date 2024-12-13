"""Contains helper functions for multiprocessing"""


# from python
from multiprocessing import Process, Pipe
from multiprocessing.connection import Connection
from typing import Callable
import logging


def start_processes(
    num_processes: int,
    worker_task_method: Callable,
    extra_data: object,
    use_batches: bool = True,
) -> tuple[list[Process], list[Connection]]:
    """Start a number of processes and return a list of processes and connections

    Args:
        num_processes (int): Number of processes to create
        worker_task_method (callable): The method that is called inside of
        worker_template
        extra_data (object): extra data that is passed to the workers. This should not
        change while the workers are active

    Returns:
        tuple[list[Process], list[Connection]]: List of processes and connections
    """
    processes: list[Process] = []
    connections: list[Connection] = []
    for process_id in range(num_processes):
        main_connection, worker_connection = Pipe()
        connections.append(main_connection)

        process = Process(
            target=worker_method,
            args=(
                process_id,
                worker_connection,
                worker_task_method,
                extra_data,
                use_batches,
            ),
        )
        processes.append(process)
        process.start()

    return processes, connections


def worker_method(
    process_id: int,
    worker_connection: Connection,
    worker_task_method: Callable,
    extra_data: object,
    use_batches=True,
) -> None:
    """Method for multiprocessing work. This function creates a loop waiting for data
    to come in on a connection, sends that data to a method, and returns the output of
    that method to the main connection

    If use_batches is set to true, data is expected to be either None or a list where the
    first element is an integer indicating the size of the task batch, and the rest of
    the elements are objects that are sent to the worker_task_method.
    If use_batches is set to false, data can be any object that is sent directly to the
    worker_task_method
    In both cases if data is none, this worker terminates.

    Args:
        process_id (int): ID of this worker
        worker_connection (Connection): connection to the main thread
        worker_task_method (callable): Method to be executed. will be passed task data
        and extra data as arguments. Must accept two parameters: task and extra_data.
        Task may change, but extra_data does not
        extra_data (any): extra data that is passed to the worker. This should not
        change while the worker is active
        use_batches (bool): True if this method uses batches, false if not
    """

    logging.debug(f"Starting worker {process_id}")

    worker_connection.send(None)

    while True:
        task_data = worker_connection.recv()

        if task_data is None:
            return

        if use_batches:
            output = [worker_task_method(task, extra_data) for task in task_data]

            worker_connection.send(output)
            continue

        output = worker_task_method(task_data, extra_data)

        worker_connection.send(output)
