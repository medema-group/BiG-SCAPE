"""Contains methods to press an input .hmm file into more optimal formats for hmmscan and hmmalign
"""

# from python
from multiprocessing.connection import Connection
from pathlib import Path
import logging

# from dependencies
from pyhmmer.plan7 import HMMFile, Pipeline, Profile
from pyhmmer.hmmer import hmmpress
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock

# from other modules
from src.genbank import CDS
from src.multithreading import Worker


class HMMer:
    pipeline: Pipeline = None
    profiles: list[Profile]

    @staticmethod
    def press(hmm_path: Path) -> None:
        """Presses hmm files for optimized files for further analysis

        Args:
            hmm_path (Path): Path to the .hmm file containing domain models

        Raises:
            ValueError: Raised if the given HMM file cannot be read by pyhmmer
        """
        with HMMFile(hmm_path) as hmm_file:
            if hmm_file.is_pressed():
                logging.info("PFAM Hmm file was already pressed.")
                return

            logging.info("Pressing Pfam HMM file...")
            hmmpress(hmm_file, hmm_path)

    @staticmethod
    def init(hmm_path: Path) -> None:
        with HMMFile(hmm_path) as hmm_file:
            HMMer.profiles = []
            HMMer.profiles.extend(list(hmm_file.optimized_profiles()))
            logging.info("Found %d hmm profiles", len(HMMer.profiles))

        # this pipeline is used to perform
        HMMer.pipeline = Pipeline(
            Alphabet.amino(), Z=len(HMMer.profiles), bit_cutoffs="trusted"
        )

    @staticmethod
    def scan_worker_method(id: int, connection: Connection) -> None:
        """Worker method for workers performing HMMScan on a sequence

        Args:
            id (int): the id of this worker
            connection (Connection): the connection through which to send and receive
            data
        """
        try:
            # we will need this many times later on
            alphabet = Alphabet.amino()
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
                    # we need the relevant data to do a scan: the id of the CDS and
                    # the protein sequence. this is what is passed into the task data
                    task_data = message[1:]
                    cds_id, aa_seq = task_data
                    # now we need to convert the sequence into something pyhmmer expects
                    text_sequence = TextSequence(sequence=aa_seq)
                    sequence_block = TextSequenceBlock([text_sequence]).digitize(
                        alphabet
                    )

                    for profile in HMMer.profiles:
                        search_result = HMMer.pipeline.search_hmm(
                            profile, sequence_block
                        )

                        for hit in search_result:
                            if not hit.reported:
                                continue

                            for domain in hit.domains.reported:
                                response = (
                                    Worker.RESULT,
                                    cds_id,
                                    domain.alignment.hmm_name,
                                    domain.score,
                                )
                                connection.send(response)

        except Exception:
            connection.send((Worker.ERROR, id))

    @staticmethod
    def scan(sequences: list[CDS]) -> None:
        """Performs hmmscan on all CDS in the database using the provided pyhmmer
        pipeline. This uses a workerpool to divide work over the CPU cores

        Args:
            pipeline (pyhmmer.plan7.pipeline): pyhmmer pipeline object created in
            init_pipeline
        """
        return
