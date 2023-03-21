"""Contains tests to test the worker/multithreading implementation in BiG-SCAPE"""

# from python
from multiprocessing.connection import Connection
from typing import Callable
from unittest import TestCase
from multiprocessing import cpu_count
from random import randint
from math import sqrt

# from other modules
from src.multithreading import WorkerPool, Worker
from src.errors import WorkerPoolSetupError

TEST_INPUT_NUM = 1000


class TestWorkerPool(TestCase):
    """Contains tests to cover worker pool functionality"""

    @staticmethod
    def kill_processes(workers: list[Worker]):
        for worker in workers:
            if worker.process is None:
                continue

            if not worker.process.is_alive():
                continue

            worker.process.kill()

    def cleanUp(self):
        """Ensures that all processes are killed at the end of testing"""
        for worker_pool in self.worker_pools:
            if len(worker_pool.workers) == 0:
                continue

            TestWorkerPool.kill_processes(worker_pool.workers)

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

    def test_start_pool(self):
        """Tests creation and start of a worker pool"""

        worker_pool = WorkerPool()
        self.worker_pools.append(worker_pool)

        def worker_function(connection: Connection):
            while True:
                command = connection.recv()
                if command[0] == Worker.COMM_STOP:
                    return

        worker_pool.worker_function = worker_function

        worker_pool.start()

        expected_worker_count = cpu_count()

        # in order to get the actual revelant processes, we just go through the workers
        # in the pool and check if they all have an active process
        worker_processes = [worker.process for worker in worker_pool.workers]

        # filter the processes to those that are alive
        alive_processes = filter(lambda process: process.is_alive(), worker_processes)
        actual_worker_count = len(list(alive_processes))

        self.assertEqual(expected_worker_count, actual_worker_count)

    def test_start_pool_no_function(self):
        """Tests whether the test pool creation throws a WorkerPoolSetupError when it is
        started without being assigned a function
        """

        worker_pool = WorkerPool()
        self.worker_pools.append(worker_pool)

        self.assertRaises(WorkerPoolSetupError, worker_pool.start)

    def test_stop_pool(self):
        """Tests whether the pool is correctly cleaned up when stop() is called"""
        worker_pool = WorkerPool()
        self.worker_pools.append(worker_pool)

        def worker_function(connection: Connection):
            while True:
                command = connection.recv()
                if command[0] == Worker.COMM_STOP:
                    return

        worker_pool.worker_function = worker_function

        worker_pool.start()

        worker_pool.stop()

        worker_pool.join()

        expected_worker_count = 0
        # in order to get the actual revelant processes, we just go through the workers
        # in the pool and check if they all have an active process
        worker_processes = [worker.process for worker in worker_pool.workers]

        # filter the processes to those that are alive
        alive_processes = filter(lambda process: process.is_alive(), worker_processes)
        actual_worker_count = len(list(alive_processes))

        self.assertEqual(expected_worker_count, actual_worker_count)

    def test_queue_task(self):
        """Tests whether a task can be successfully queued on a worker pool"""
        pass

    def test_set_callback(self):
        """Tests whether a callback can be successfully set on a worker pool"""

    def test_set_callback_already_set(self):
        """Tests whether set_callback raises a WorkerPoolSetupError when a callback is
        already set and a user tries to set another one
        """


class TestWorker(TestCase):
    """Contains tests to cover worker functionality"""

    @staticmethod
    def kill_processes(workers: list[Worker]):
        for worker in workers:
            if worker.process is None:
                continue

            if not worker.process.is_alive():
                continue

            worker.process.kill()

    def cleanUp(self):
        """Ensures that all processes are killed at the end of testing"""
        TestWorker.kill_processes(self.workers)

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        # keep track of workers to clean up later
        self.workers: list[Worker] = []

        self.addCleanup(self.cleanUp)

    def test_start_worker(self):
        """Tests whether starting a worker spawns a subprocess"""

        # we are not using the worker pool here beyond supplying a function
        worker_pool = WorkerPool(0)

        def worker_function(connection: Connection):
            while True:
                command = connection.recv()
                if command[0] == Worker.COMM_STOP:
                    return

        worker_pool.worker_function = worker_function

        worker = Worker(0, worker_pool)
        self.workers.append(worker)

        worker.start()

        self.assertTrue(worker.process.is_alive())

    def test_stop_worker(self):
        """Tests whether calling stop on worker actually stops the worker"""

        # we are not using the worker pool here beyond supplying a function
        worker_pool = WorkerPool(0)

        def worker_function(connection: Connection):
            while True:
                command = connection.recv()
                if command[0] == Worker.COMM_STOP:
                    return

        worker_pool.worker_function = worker_function

        worker = Worker(0, worker_pool)
        self.workers.append(worker)

        worker.start()

        worker.stop()

        worker.join()

        self.assertFalse(worker.process.is_alive())
