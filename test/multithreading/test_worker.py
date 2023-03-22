"""Contains tests to test the worker/multithreading implementation in BiG-SCAPE"""

# from python
from multiprocessing.connection import Connection, Pipe
from unittest import TestCase
from multiprocessing import cpu_count

# from other modules
from src.multithreading import WorkerPool, Worker
from src.errors import WorkerPoolSetupError


class TestWorkerPool(TestCase):
    """Contains tests to cover worker pool functionality"""

    @staticmethod
    def worker_stop_method(id: int, connection: Connection):
        """Worker function that immediately stops after receiving any command"""
        while True:
            command = connection.recv()
            if command[0] == Worker.STOP:
                return

    def worker_echo(id, connection: Connection):
        """Proper implementation of a worker function which echos any task sent
        to it back as a result
        """

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
            TestWorkerPool.kill_processes(worker_pool.workers)

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        # add clean up method
        self.addCleanup(self.cleanUp)

        # we will collect the worker pools here so that we can make sure all processes
        # are closed off
        self.worker_pools: list[WorkerPool] = []

    def test_worker_stop_method(self):
        """Tests the simple worker method built into this test class"""
        test_connection, process_connection = Pipe(True)

        command = (Worker.STOP, 0)
        test_connection.send(command)

        TestWorkerPool.worker_stop_method(0, process_connection)

        self.assertTrue(test_connection.writable)

    def test_start_pool(self):
        """Tests creation and start of a worker pool"""

        worker_pool = WorkerPool()
        self.worker_pools.append(worker_pool)

        worker_pool.worker_function = TestWorkerPool.worker_stop_method

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

        worker_pool.worker_function = TestWorkerPool.worker_stop_method

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


class TestWorker(TestCase):
    """Contains tests to cover worker functionality"""

    @staticmethod
    def kill_processes(workers: list[Worker]):
        for worker in workers:
            if worker.process is None or not worker.process.is_alive():
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

        worker_pool.worker_function = TestWorkerPool.worker_stop_method

        worker = Worker(0, worker_pool)
        self.workers.append(worker)

        pool_connection, worker_connection = Pipe(True)
        worker.start(pool_connection, worker_connection)

        self.assertTrue(worker.process.is_alive())

    def test_stop_worker(self):
        """Tests whether calling stop on worker actually stops the worker"""

        # we are not using the worker pool here beyond supplying a function
        worker_pool = WorkerPool(0)

        worker_pool.worker_function = TestWorkerPool.worker_stop_method

        worker = Worker(0, worker_pool)
        self.workers.append(worker)

        # workers need connections
        pool_connection, worker_connection = Pipe(True)
        worker.start(pool_connection, worker_connection)

        worker.stop()

        worker.join()

        self.assertFalse(worker.process.is_alive())

    def test_kill_worker(self):
        """Tests whether calling stop on worker actually stops the worker"""

        # we are not using the worker pool here beyond supplying a function
        worker_pool = WorkerPool(0)

        worker_pool.worker_function = TestWorkerPool.worker_stop_method

        worker = Worker(0, worker_pool)
        self.workers.append(worker)

        # workers need connections
        pool_connection, worker_connection = Pipe(True)
        worker.start(pool_connection, worker_connection)

        worker.stop(True)

        worker.join()

        self.assertFalse(worker.process.is_alive())
