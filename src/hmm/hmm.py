"""Contains methods to press an input .hmm file into more optimal formats for hmmscan and hmmalign
"""

# from python
import logging
from typing import Iterator, cast
from math import ceil
from pathlib import Path
from multiprocessing import Pipe, Process, cpu_count
from multiprocessing.connection import Connection, wait

# from dependencies
from pyhmmer.plan7 import (
    HMMFile,
    Pipeline,
    OptimizedProfile,
)
from pyhmmer.hmmer import hmmpress, hmmsearch
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock


# from other modules
from src.genbank import CDS

# from this module
from src.hmm.hsp import HSP


class HMMer:
    pipeline: Pipeline = None
    profiles: list[OptimizedProfile]
    alphabet = Alphabet.amino()

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
        logging.info("Reading HMM Profiles")
        HMMer.profiles = list()
        with HMMFile(hmm_path) as hmm_file:
            HMMer.profiles.extend(list(hmm_file.optimized_profiles()))
            logging.info("Found %d HMM profiles", len(HMMer.profiles))

        # this pipeline is used to perform
        HMMer.pipeline = Pipeline(
            Alphabet.amino(), Z=len(HMMer.profiles), bit_cutoffs="trusted"
        )

    @staticmethod
    def search(genes: list[CDS]) -> Iterator[HSP]:
        """Performs hmmsearch on a list of CDS.

        This is the fastest method of getting scan results. Use this if you know you
        won't run into memory issues.

        Args:
            genes (list[CDS]): list of genes to generate HSPs for

        Yields:
            Iterator[HSP]: an iterator generating HSPs that are included
        """

        logging.info("Performing search on %d genes", len(genes))

        sequences = []
        for idx, gene in enumerate(genes):
            sequences.append(TextSequence(name=str(idx).encode(), sequence=gene.aa_seq))

        ds = TextSequenceBlock(sequences).digitize(HMMer.alphabet)

        for top_hits in hmmsearch(HMMer.profiles, ds, bit_cutoffs="trusted"):
            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    if domain.score < 0:
                        continue

                    cds_idx = int(hit.name.decode())
                    accession = domain.alignment.hmm_accession.decode()
                    score = domain.score
                    yield HSP(genes[cds_idx], accession, score)

    @staticmethod
    def scan(genes: list[CDS], batch_size=None) -> Iterator[HSP]:
        """Runs hmmscan using pyhmmer in several subprocesses. Passes batches of input
        genes to the subprocesses in an attempt to optimize memory usage

        Args:
            genes (list[CDS]): a list of CDS objects to generate HSPs for
            batch_size (int, optional): The size of the data batch lists which is passed
            to the worker processes. Defaults to None.

        Yields:
            Iterator[HSP]: An iterator generating hsps from the list of CDS passed into
            this method
        """

        processes = []
        connections = []

        if batch_size is None:
            batch_size = ceil(len(genes) / cpu_count())
            logging.info("Using automatic batch size %d", batch_size)
        else:
            logging.info("Using manual batch size %d", batch_size)

        task_iter = task_generator(genes, batch_size)

        # first we need to create the processes and the connections
        for process_id in range(cpu_count()):
            # connection first. this is the communication between the main process
            # and the worker process
            main_connection, worker_connection = Pipe(True)

            # we need to keep track of connections in the coming code for sending and
            # receiving
            connections.append(main_connection)

            # now we can create a process
            # we need to tell it what method to execute, and give it some parameters
            # one of the arguments is the other end of the connection we made before
            worker_process = Process(
                target=HMMer.hmmscan_process, args=(process_id, worker_connection)
            )
            processes.append(worker_process)

            # start the process
            worker_process.start()

        # the trouble with the above method of doing multiprocessing is that the
        # connections are a big hassle. it's easy to end up in a situation where some
        # process or the main thread gets locked, e.g. when both ends are trying to read
        # or both ends are trying to write
        # in order to avoid any trouble we will try to make it so that the workers will
        # always be in a waiting state by trying to send something to the main process
        # the main process will receive this and send back either a set of tasks or
        # an order to shut down
        # this way we can write the main process as follows:

        tasks_done = 0

        while len(connections) > 0:
            # get the connections (workers) which are reable (trying to send us stuff)
            available_connections = wait(connections)

            # go through all of these
            for connection in available_connections:
                connection = cast(Connection, connection)
                # receive the data and store for now. we will handle it later
                output_data = connection.recv()

                # see if there is any data to get
                input_data = next(task_iter, None)

                # if there is no data, input_data will be set to None
                # this tells the worker to shut down
                # if there is, it's a list that can be sent to the worker directly
                connection.send(input_data)

                # we do need to clean up the connection if it is going to close, though
                if input_data is None:
                    connection.close()
                    connections.remove(connection)

                # now lets handle the output data
                # it will be either None or a list
                # the fist element of the list will be a number indicating the amount of
                # tasks this output came from
                # any other elements are a tuple of format:
                # (cds_id, score, domain_accession)
                # we only care about if a list is returned. otherwise we just move on
                if output_data is not None:
                    src_task_count: int = output_data[0]
                    task_outputs: list[tuple[int, str, float]] = output_data[1:]
                    for task_output in task_outputs:
                        yield task_output_to_hsp(task_output, genes)

                    tasks_done += src_task_count

    @staticmethod
    def profile_hmmsearch(
        sequences: list[TextSequence],
    ) -> list[tuple[int, str, float]]:
        """Performs search_hmm on the given list of sequences for each profile in
        HMMer.profiles

        Args:
            sequences (list[TextSequence]): list of input sequences

        Returns:
            list[tiple[int, str, float]]: result object to pass back to the main thread
        """
        text_seq_block = TextSequenceBlock(sequences).digitize(HMMer.alphabet)

        outputs = []
        # perform hmmsearch for each profile
        for profile in HMMer.profiles:
            hmmscan_results = HMMer.pipeline.search_hmm(profile, text_seq_block)
            for top_hit in hmmscan_results.reported:
                for domain in top_hit.domains:
                    if not domain.reported:
                        continue

                    cds_idx = int(top_hit.name.decode())
                    accession = domain.alignment.hmm_accession.decode()
                    score = domain.score
                    outputs.append((cds_idx, accession, score))

        return outputs

    @staticmethod
    def hmmscan_process(process_id, connection: Connection):  # pragma: no cover
        """Process for hmmscan workers

        Args:
            process_id (int): the id of this process
            connection (multiprocessing.connection.Connection): Worker connection of the
                pipe connecting to the main process
        """
        # start by waiting for the main thread to be ready for us
        connection.send(None)
        while True:
            # now we can get down to business
            # get a batch of tasks
            tasks = connection.recv()

            if tasks is None:
                break

            num_tasks = len(tasks)

            # get sequences
            sequences = []
            for task in tasks:
                cds_id, aa_seq = task
                sequences.append(
                    TextSequence(name=str(cds_id).encode(), sequence=aa_seq)
                )

            task_output = HMMer.profile_hmmsearch(sequences)

            connection.send([num_tasks] + task_output)


def cds_to_input_task(cds_list: list[CDS]) -> Iterator[tuple[int, str]]:
    """Returns an iterator which yields input tasks for each cds in the given list

    A CDS without an amino acid sequence will be passed as an empty sequence string

    Args:
        cds_list (list[CDS]): list of CDS to generate tasks for

    Yields:
        Iterator[tuple[int, str]]: An iterator that generates tasks, which are tuples of
        the format (index: int, aa_seq: str)
    """
    for idx, cds in enumerate(cds_list):
        yield (idx, cds.aa_seq or "")


def task_output_to_hsp(task_output: tuple, cds_list: list[CDS]):
    """Returns an HSP for a given output task"""
    cds_id, domain, score = task_output
    return HSP(cds_list[cds_id], domain, score)


def task_generator(genes: list[CDS], batch_size) -> Iterator[list]:
    """Returns an iterator which yields subsets of a given list, split by batch size

    Args:
        full_list (list): a larger list
        batch_size (int): the maximum size of a list chunk
    """
    batches, remainder = divmod(len(genes), batch_size)
    for i in range(batches + 1):
        start = i * batch_size
        stop = start + batch_size

        if stop > len(genes):
            stop = start + remainder
        gene_batch = genes[start:stop]

        tasks = []
        for cds_idx, gene in enumerate(gene_batch):
            tasks.append((cds_idx, gene.aa_seq))

        yield tasks


def hsp_overlap_filter(hsp_list: list[HSP], overlap_cutoff=0.1) -> list[HSP]:
    """Filters overlapping HSP

    Args:
        hsp_list (list[HSP]): a list of high scoring protein hits
        overlap_cutoff: The maximum percentage sequence length domains can overlap
            without being discarded

    Returns:
        list[HSP]: a filtered list of high scoring protein hits
    """
    # for this the hsps have to be sorted by any of the start coordinates
    # this is for the rare case where two hsps have the same bit score.
    # to replicate BiG-SCAPE output, the hsp with a lower start coordinate
    # will end up being chosen
    sorted_hsps = sorted(hsp_list, key=lambda hsp: hsp.cds.nt_start)

    new_list = []
    for i in range(len(sorted_hsps) - 1):
        for j in range(i + 1, len(sorted_hsps)):
            hsp_a: HSP = sorted_hsps[i]
            hsp_b: HSP = sorted_hsps[j]

            # check if we are the same CDS
            if hsp_a == hsp_b:
                # check if there is overlap between the domains
                if not has_overlap(hsp_a, hsp_b):
                    continue

                overlapping_aminoacids = len_overlap(hsp_a, hsp_b)
                overlap_perc_loc1 = overlapping_aminoacids / (
                    hsp_a.cds.nt_stop - hsp_a.cds.nt_start
                )
                overlap_perc_loc2 = overlapping_aminoacids / (
                    hsp_b.cds.nt_stop - hsp_b.cds.nt_start
                )
                # check if the amount of overlap is significant
                if (
                    overlap_perc_loc1 <= overlap_cutoff
                    and overlap_perc_loc2 <= overlap_cutoff
                ):
                    continue
                # rounding with 1 decimal to mirror domtable files
                hsp_a_score = round(float(hsp_a.score), 1)
                hsp_b_score = round(float(hsp_b.score), 1)
                if hsp_a_score >= hsp_b_score:  # see which has a better score
                    new_list.append(hsp_a)
                elif hsp_a_score < hsp_b_score:
                    new_list.append(hsp_b)

    return new_list


def has_overlap(hsp_a: HSP, hsp_b: HSP):
    """Return True if there is overlap between two regions"""
    # a left of b
    if hsp_a.cds.nt_stop < hsp_b.cds.nt_start:
        return False
    # b left of a
    if hsp_b.cds.nt_stop < hsp_a.cds.nt_start:
        return False

    # all other cases should have overlap
    return True


def len_overlap(hsp_a: HSP, hsp_b: HSP):
    """Returns the length of an overlapping sequence"""

    if hsp_a.cds.nt_start < hsp_b.cds.nt_start:
        left = hsp_b.cds.nt_start
    else:
        left = hsp_a.cds.nt_start

    if hsp_a.cds.nt_stop > hsp_b.cds.nt_stop:
        right = hsp_b.cds.nt_stop
    else:
        right = hsp_a.cds.nt_stop

    # limit to > 0
    return max(0, right - left)
