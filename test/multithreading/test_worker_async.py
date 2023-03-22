"""Contains tests to test the worker/multithreading implementation in BiG-SCAPE"""

# from python
import logging
from multiprocessing.connection import Connection, Pipe
from typing import Callable
from unittest import TestCase
from random import randint
from math import sqrt

# from other modules
from src.multithreading import WorkerPool, Worker

TEST_INPUT_NUM = 1000


class TestWorkerPoolAsync(TestCase):
    """Contains tests to cover asynchronous parts of test workers"""

    def worker_echo_method(id, connection: Connection):
        try:
            while True:
                message = connection.recv()
                command_type = message[0]

                logging.debug("W%d: %d", id, message[0])

                if command_type == Worker.START:
                    connection.send(message)
                    continue

                if command_type == Worker.STOP:
                    connection.send(message)
                    break

                if command_type == Worker.TASK:
                    response = (Worker.RESULT,) + message[1:]
                    connection.send(response)
        except Exception:
            connection.send((Worker.ERROR, id))

    @staticmethod
    def kill_processes(workers: list[Worker]):
        for worker in workers:
            if worker.process is None or not worker.process.is_alive():
                continue

            worker.process.kill()

    def cleanUp(self):
        """Ensures that all processes are killed at the end of testing"""
        for worker_pool in self.worker_pools:
            TestWorkerPoolAsync.kill_processes(worker_pool.workers)

    def setUp(self) -> None:
        """Contains setup necessary for testing"""
        # generate a set of test inputs to use later as an expected and actual output
        for i in range(TEST_INPUT_NUM):
            number = randint(0, 100)
            answer = self.test_function(number)

            self.test_inputs.append((i, number))
            self.test_outputs.append((i, answer))

        return super().setUp()

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        # add clean up method
        self.addCleanup(self.cleanUp)

        # the method we will pass to our workers later on
        def test_function(number):
            return sqrt(number)

        self.test_function: Callable = test_function

        # the expected inputs/outputs
        self.test_inputs: list[tuple] = []
        self.test_outputs: list[tuple] = []

        # we will collect the worker pools here so that we can make sure all processes
        # are closed off
        self.worker_pools: list[WorkerPool] = []

    def test_worker_echo_method(self):
        """Tests the worker echo method that is part of this test class"""
        test_connection, process_connection = Pipe(True)

        expected_responses = []

        start_command = (Worker.START, 0)
        test_connection.send(start_command)
        expected_responses.append(start_command)

        # task command
        task_command = (Worker.TASK, "test")
        test_connection.send(task_command)
        expected_result = (Worker.RESULT, "test")
        expected_responses.append(expected_result)

        # stop command
        stop_command = (Worker.STOP, 0)
        test_connection.send(stop_command)
        expected_responses.append(stop_command)

        TestWorkerPoolAsync.worker_echo_method(0, process_connection)

        actual_responses = []
        # start
        actual_responses.append(test_connection.recv())

        # result
        actual_responses.append(test_connection.recv())

        # stop
        actual_responses.append(test_connection.recv())

        self.assertEqual(expected_responses, actual_responses)

    def test_worker_task(self):
        """Tests the execution of a single task by a single worker"""
        worker_pool = WorkerPool(1)
        self.worker_pools.append(worker_pool)

        # we will try to send something simple and expect the exact same thing back

        worker_pool.worker_function = TestWorkerPoolAsync.worker_echo_method

        worker_pool.start()

        expected_data = "test"

        task = (expected_data,)
        worker_pool.send_task_one(task)

        # skip START response
        worker_pool.read_one()

        # task result
        response = worker_pool.read_one()

        actual_data = response[1]

        worker_pool.stop()

        worker_pool.join()

        self.assertEqual(expected_data, actual_data)
