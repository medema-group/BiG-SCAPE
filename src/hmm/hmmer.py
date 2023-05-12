"""Contains methods to press an input .hmm file into more optimal formats for hmmscan and hmmalign
"""

# from python
from __future__ import annotations
import logging
import string
from math import ceil
from typing import Callable, Iterator, Optional, cast, TYPE_CHECKING
from pathlib import Path
from multiprocessing import Pipe, Process, cpu_count
from multiprocessing.connection import Connection, wait

# from dependencies
from pyhmmer.plan7 import (
    HMMFile,
    Pipeline,
    OptimizedProfile,
)
from pyhmmer.hmmer import hmmpress, hmmsearch, hmmalign
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock


# from this module
from src.hmm.hsp import HSP, HSPAlignment

# circular dependencies
if TYPE_CHECKING:  # pragma: no cover
    from src.genbank import CDS  # imported in genbank.GBK


class HMMer:
    pipeline: Pipeline = None
    profiles: list[OptimizedProfile]
    profile_index: dict[str, int]
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
    def init(hmm_path: Path, optimized=True) -> None:
        """Load the HMMer profiles and create an index for the profiles for easy lookup
        later. Also creates a pipeline for usage in certain workflows

        Args:
            hmm_path (Path): Path to the hmm profile. Should contain hmm[f,i,m,p] if
            optimized is set to True
            optimized (bool, optional): Whether to use the optimized profiles if they
            exist. Defaults to True.
        """
        logging.info("Reading HMM Profiles")
        HMMer.profiles = list()
        with HMMFile(hmm_path) as hmm_file:
            if optimized:
                logging.debug("Loading optimized profiles")
                profile_list = list(hmm_file.optimized_profiles())
            else:
                logging.debug("Loading unoptimized profiles")
                profile_list = list(hmm_file)

            HMMer.profiles.extend(profile_list)
            logging.info("Found %d HMM profiles", len(HMMer.profiles))

        # create an index so we can quickly get profiles by accession
        HMMer.profile_index = gen_profile_index(HMMer.profiles)

        # this pipeline is used to perform hmmscan and hmmalign in more memory efficient
        # use cases. This pipeline functions as a wrapper for hmmsearch etc. and can be
        # used to pass arguments that you can otherwise pass to these functions
        HMMer.pipeline = Pipeline(
            Alphabet.amino(), Z=len(HMMer.profiles), bit_cutoffs="trusted"
        )

    @staticmethod
    def unload():
        """Unload the variables that were set up during init(). Do this before running
        init() again for hmmalign
        """
        logging.info("Unloading profiles")
        HMMer.profiles = list()
        HMMer.profile_index = {}
        HMMer.pipeline = None

    @staticmethod
    def get_profile(accession: str) -> OptimizedProfile:
        """Returns the pressed/optimized profile with the specified accession

        Args:
            accession (str): the accession identifying the profile of interest

        Returns:
            OptimizedProfile: An optimized/pressed HMM profile
        """
        if accession not in HMMer.profile_index:
            raise KeyError(
                "KeyError in retrieving profile. Was a different .hmm used "
                "during alignment, or did the .hmm file change during "
                "execution?"
            )

        profile_idx = HMMer.profile_index[accession]
        return HMMer.profiles[profile_idx]

    @staticmethod
    def hmmsearch_simple(
        cds_list: list[CDS], domain_overlap_cutoff=0.1, cores=cpu_count()
    ):
        """Performs hmmsearch on a list of CDS and appends detected HSPs to the CDS
        objects

        HSP objects contain domain and score information, and a link to the related CDS

        Args:
            cds_list (list[CDS]): list of CDS to generate HSPs for
            domain_overlap_cutoff (float): cutoff to use for domain overlap filtering.
            Defaults to 0.1
        """

        logging.info("Performing simple hmmsearch on %d genes", len(cds_list))

        sequences = []
        for idx, cds in enumerate(cds_list):
            # name is actually the list index of the original CDS
            sequences.append(TextSequence(name=str(idx).encode(), sequence=cds.aa_seq))

        ds = TextSequenceBlock(sequences).digitize(HMMer.alphabet)

        for top_hits in hmmsearch(
            HMMer.profiles, ds, bit_cutoffs="trusted", cpus=cores
        ):
            for hit in top_hits:
                if not hit.included:
                    continue
                for domain in hit.domains:
                    if not domain.included:
                        continue
                    if domain.score < 0:
                        continue

                    # name is actually the list index of the original CDS
                    cds_idx = int(hit.name.decode())
                    accession = domain.alignment.hmm_accession.decode()
                    score = domain.score
                    env_start = domain.env_from
                    env_stop = domain.env_to
                    relevant_cds = cds_list[cds_idx]
                    hsp = HSP(relevant_cds, accession, score, env_start, env_stop)
                    relevant_cds.add_hsp_overlap_filter(hsp, domain_overlap_cutoff)

    @staticmethod
    def calc_batch_size(num_genes: int, cores=cpu_count()) -> int:
        """calculate batch size by splitting total tasks/genes between available CPUs

        Args:
            num_genes (int): number of genes that are to be processed

        Returns:
            int: the gene batch size for worker processes
        """
        return ceil(num_genes / cores)

    @staticmethod
    def hmmsearch_multiprocess(
        cds_list: list[CDS],
        domain_overlap_cutoff: float = 0.1,
        cores=cpu_count(),
        batch_size: Optional[int] = None,
        callback: Optional[Callable] = None,
    ) -> None:
        """Runs hmmscan using pyhmmer in several subprocesses. Passes batches of input
        cds list to the subprocesses in an attempt to optimize memory usage

        Args:
            cds_list (list[CDS]): a list of CDS objects to generate HSPs for
            domain_overlap_cutoff (float): cutoff to use for domain overlap filtering.
            Defaults to 0.1
            cores (int): number of processes to use. Defaults to number of logical cores
            batch_size (int, optional): The size of the data batch lists which is passed
            to the worker processes. Defaults to None.
            callback (Callable, optional): A callback function with a num_tasks argument
            that is called whenever a set of tasks is done

        """
        logging.info("Performing distributed hmmsearch on %d genes", len(cds_list))
        processes: list[Process] = []
        connections: list[Connection] = []

        # check for set batch size
        if batch_size is None:
            batch_size = HMMer.calc_batch_size(len(cds_list), cores)
            logging.info("Using automatic batch size %d", batch_size)

        task_iter = task_generator(cds_list, batch_size)

        # it doesn't make sense to spawn more processes than there are AA sequences to
        # search, so only spawn between 0 and cpu count processes
        process_count = min(len(cds_list), cores)

        # first we need to create the worker processes and the connections
        for process_id in range(process_count):
            # connection first. this is the communication between the main process
            # and the worker process
            # note that here worker_connection is a connection the worker will use
            # and main_connection is the connection the main process will use
            main_connection, worker_connection = Pipe(True)

            # we need to keep track of connections in the coming code for sending and
            # receiving to and from worker processes
            connections.append(main_connection)

            # now we can create a worker process
            # we need to tell it what method to execute, and give it some parameters.
            # one of the arguments is the other end of the connection we made before
            worker_process = Process(
                target=HMMer.hmmsearch_process, args=(process_id, worker_connection)
            )
            processes.append(worker_process)

            # start the worker process
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

        # now that we have our worker processes and connections, we can start sending
        # data to the workers so they can actually start doing stuff
        while len(connections) > 0:
            # get the connections which are readable. these are the workers which
            # have sent us something and are therefore able to receive something
            available_connections = wait(connections)

            # go through all of the available main connections
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

                    # list of tuples that are outputs for the tasks sent to this process
                    # tuples are in the format (cds_idx, domain_accession, score)
                    task_outputs: list[tuple[int, str, float]] = output_data[1:]
                    for task_output in task_outputs:
                        process_task_output(
                            task_output, cds_list, domain_overlap_cutoff
                        )

                    tasks_done += src_task_count
                    if callback is not None:
                        callback(tasks_done)

        # just to make sure, kill any remaining processes
        for process in processes:
            process.kill()

    @staticmethod
    def profile_hmmsearch(
        sequences: list[TextSequence],
    ) -> list[tuple[int, str, float, int, int]]:
        """Performs search_hmm on the given list of sequences for each profile in
        HMMer.profiles

        return value is a list of tuples with the structure:
        (cds_idx, domain_accession, score)
        These values are passed back to the main process

        Args:
            sequences (list[TextSequence]): list of input sequences

        Returns:
            list[tuple[int, str, float]]: result object to pass back to the main process
        """
        text_seq_block = TextSequenceBlock(sequences).digitize(HMMer.alphabet)

        outputs = []
        # perform hmmsearch for each profile
        for profile in HMMer.profiles:
            hmmsearch_results = HMMer.pipeline.search_hmm(profile, text_seq_block)
            for top_hit in hmmsearch_results.reported:
                for domain in top_hit.domains:
                    if not domain.reported:
                        continue

                    # name is actually the list index of the original CDS
                    cds_idx = int(top_hit.name.decode())
                    accession = domain.alignment.hmm_accession.decode()
                    score = domain.score
                    env_start = domain.env_from
                    env_stop = domain.env_to
                    outputs.append((cds_idx, accession, score, env_start, env_stop))

        return outputs

    @staticmethod
    def hmmsearch_process(
        process_id, worker_connection: Connection
    ):  # pragma: no cover
        """Process for hmmsearch workers

        Args:
            process_id (int): the id of this process
            connection (multiprocessing.connection.Connection): Worker connection of the
                pipe connecting to the main process
        """
        # start by waiting for the main thread to be ready for us
        logging.debug("hmmsearch process with id %d started", process_id)
        worker_connection.send(None)
        while True:
            # now we can get down to business
            # get a batch of tasks
            tasks = worker_connection.recv()

            # we established that receiving None means we should stop the process
            if tasks is None:
                logging.debug("hmmsearch process with id %d stopping", process_id)
                break

            # we will want to return the number of tasks that were processed later
            num_tasks = len(tasks)

            # pymmer requires sequences to be in a certain format:
            sequences = []
            for task in tasks:
                cds_id, aa_seq = task
                # name is actually the list index of the original CDS
                sequences.append(
                    TextSequence(name=str(cds_id).encode(), sequence=aa_seq)
                )

            task_output = HMMer.profile_hmmsearch(sequences)

            # send this batch of task results back as a list
            # first element is the number of tasks that were processed in this batch
            worker_connection.send([num_tasks] + task_output)

    @staticmethod
    def align_simple(hsps: list[HSP]) -> Iterator[HSPAlignment]:
        """Aligns a list of HSPs per domain.

        The iterator that is returned will yield a list of alignments per domain

        This performs no multiprocessing, but align is usually pretty fast already

        Args:
            hsps (list[HSP]): List of HSPs to align

        Yields:
            Iterator[HSPAlignment]: A generator of HSPAlignments
        """

        # there are going to be cases where certain domains are not present in HSPs at
        # all. we can assemble a dictionary which contains only those domains that are
        # present and the HSPs associated with those domains
        domain_hsps: dict[str, list[HSP]] = {}

        for hsp in hsps:
            domain_accession = hsp.domain
            if domain_accession not in domain_hsps:
                domain_hsps[hsp.domain] = []

            domain_hsps[hsp.domain].append(hsp)

        # now all we have to do is go through the list of domains which we know need to
        # be aligned with only the relevant aa sequences
        for domain_accession, hsp_list in domain_hsps.items():
            domain_profile = HMMer.get_profile(domain_accession)
            sequences = []

            domain_hsp: HSP
            for hsp_idx, domain_hsp in enumerate(hsp_list):
                # we can cut the alignment down to only the range of the domain as found
                # in hmmscan to speed up the process
                start = domain_hsp.env_start
                stop = domain_hsp.env_stop
                seq_to_align = domain_hsp.cds.aa_seq[start:stop]
                sequences.append(
                    TextSequence(name=str(hsp_idx).encode(), sequence=seq_to_align)
                )
            ds_block = TextSequenceBlock(sequences).digitize(HMMer.alphabet)

            msa = hmmalign(domain_profile, ds_block)

            for idx, alignment in enumerate(msa.alignment):
                # name is actually the list index of the original HSP
                sequence_name = msa.sequences[idx].name.decode()
                source_hsp_idx = int(sequence_name)

                algn_string = process_algn_string(alignment)
                yield HSPAlignment(hsp_list[source_hsp_idx], algn_string)


# general methods
def gen_profile_index(profiles: list[OptimizedProfile]) -> dict[str, int]:
    """Generates a profile index where keys are accessions and values are list indexes

    Used in quick lookup of profiles during hmmalign

    Args:
        profiles (list[OptimizedProfile]): a list of optimized profiles

    Returns:
        dict[str, int]: An index dictionary allowing for lookup in list by accession
    """
    index = {}
    for idx, profile in enumerate(profiles):
        index[profile.accession.decode()] = idx

    return index


def cds_to_input_task(cds_list: list[CDS]) -> Iterator[tuple[int, str]]:
    """Returns an iterator which yields input tasks for each CDS in the given list

    A CDS without an amino acid sequence will be passed as an empty sequence string

    Args:
        cds_list (list[CDS]): list of CDS to generate tasks for

    Yields:
        Iterator[tuple[int, str]]: An iterator that generates tasks, which are tuples of
        the format (index: int, aa_seq: str)
    """
    for idx, cds in enumerate(cds_list):
        yield (idx, cds.aa_seq or "")


def process_task_output(
    task_output: tuple, cds_list: list[CDS], domain_overlap_cutoff
) -> None:
    """Adds an HSP from a task output to the relevant CDS

    Args:
        task_output (tuple): an output task as given by HMMer.hmmsearch_process
        cds_list (list[CDS]): _description_
    """
    cds_id, domain, score, env_start, env_stop = task_output
    relevant_cds: CDS = cds_list[cds_id]
    hsp = HSP(relevant_cds, domain, score, env_start, env_stop)
    relevant_cds.add_hsp_overlap_filter(hsp, domain_overlap_cutoff)


def task_generator(cds_list: list[CDS], batch_size) -> Iterator[list]:
    """Returns an iterator which yields subsets of the given CDS list in task format
    for subprocesses

    tasks are in the format (cds_id, aa_seq)

    Args:
        full_list (list): a larger list
        batch_size (int): the maximum size of a list chunk

    Yields:
        Iterator[list[tuple[int, str]]]: Iterator of lists of task tuples
    """
    batches, remainder = divmod(len(cds_list), batch_size)
    for i in range(batches + 1):
        start = i * batch_size
        stop = start + batch_size

        if stop > len(cds_list):
            stop = start + remainder
        cds_batch = cds_list[start:stop]

        tasks = []
        for batch_idx, cds in enumerate(cds_batch):
            cds_idx = start + batch_idx
            tasks.append((cds_idx, cds.aa_seq))

        yield tasks


def process_algn_string(algn_string: str):
    """removes any model gaps from the alignment string and returns a string with only
    cds gaps
    """
    return algn_string.translate(str.maketrans("", "", string.ascii_lowercase + "."))
