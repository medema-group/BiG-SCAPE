"""Contains tests to test the hmm scanning functionality using pyhmmer"""

# from python
from multiprocessing import Pipe
from pathlib import Path
from unittest import TestCase

# from other modules
from src.hmm import HMMer
from src.multithreading import Worker


class TestHMMScan(TestCase):
    """Contains tests to check the hmmscan functionality"""

    def test_scan_worker_method(self):
        """Tests scanning of a single sequence on a set of domain HMMs"""

        aa_seq = (
            "MQQDGTQQDRIKQSPAPLNGMSRRGFLGGAGTLALATASGLLLPGTAHAATTITTNQTGTDGMYYSFWTDGGGS"
            "VSMTLNGGGSYSTQWTNCGNFVAGKGWSTGGRRTVRYNGYFNPSGNGYGCLYGWTSNPLVEYYIVDNWGSYRPT"
            "GTYKGTVSSDGGTYDIYQTTRYNAPSVEGTKTFQQYWSVRQSKVTSGSGTITTGNHFDAWARAGMNMGQFRYYM"
            "IMATEGYQSSGSSNITVSG"
        )

        # loading this specific hmm because we know the above sequences will be matched
        # by it
        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        test_connection, scan_method_connection = Pipe(True)

        task_command = (Worker.TASK, 1, aa_seq)

        test_connection.send(task_command)

        # send a stop command as well so that the worker can close when its done

        stop_command = (Worker.STOP, 1)
        test_connection.send(stop_command)

        HMMer.scan_worker_method(0, scan_method_connection)

        response = test_connection.recv()

        expected_response_type = Worker.RESULT
        actual_response_type = response[0]

        # just doing an assert equal. if we get a result it means we got a hit
        self.assertEqual(expected_response_type, actual_response_type)

    def test_scan_worker_method_exception(self):
        """Tests if the hmmscan worker method returns an exception message if something
        goes wrong
        """

        hmm_path = Path("test/test_data/hmm/PF00457.19.hmm")

        HMMer.init(hmm_path)

        test_connection, scan_method_connection = Pipe(True)

        # sending a command with a number instead of a string will certainly cause
        # issues
        task_command = (Worker.TASK, 1, 1)

        test_connection.send(task_command)

        # send a stop command as well so that the worker can close when its done

        stop_command = (Worker.STOP, 1)
        test_connection.send(stop_command)

        HMMer.scan_worker_method(0, scan_method_connection)

        response = test_connection.recv()

        expected_response_type = Worker.ERROR
        actual_response_type = response[0]

        # just doing an assert equal. if we get a result it means we got a hit
        self.assertEqual(expected_response_type, actual_response_type)
