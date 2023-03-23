"""Contains methods to press an input .hmm file into more optimal formats for hmmscan and hmmalign
"""

# from python
from pathlib import Path
import logging
from typing import Iterator

# from dependencies
from pyhmmer.plan7 import HMMFile, Pipeline, Profile, TopHits
from pyhmmer.hmmer import hmmpress, hmmsearch
from pyhmmer.easel import Alphabet, TextSequence, TextSequenceBlock


# from other modules
from src.genbank import CDS

# from this module
from src.hmm.hsp import HSP


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
    def scan(genes: list[CDS]) -> Iterator[HSP]:
        """Performs hmmscan on all CDS in the database using pyhmmer"""

        # callback function for progress reporting
        def callback(query, profiles_done):
            """Callback function for progress reporting

            Args:
                query (pyhmmer.plan7.OptimizedProfile): The profile that was just searched
                profiles_done (int): the number of profiles searched through so far
            """
            percentage = int(profiles_done / len(HMMer.profiles))
            report = (
                len(HMMer.profiles) < 100
                or profiles_done % int(len(HMMer.profiles) / 100) == 0
            )
            if report:
                logging.info(
                    "%d/%d (%d%%)", profiles_done, len(HMMer.profiles), percentage
                )

        text_sequences = []
        for idx, gene in enumerate(genes):
            text_sequences.append(
                TextSequence(name=str(idx).encode(), sequence=gene.aa_seq)
            )

        alphabet = Alphabet.amino()
        sequence_block = TextSequenceBlock(text_sequences).digitize(alphabet)

        hmmsearch_iter = hmmsearch(HMMer.profiles, sequence_block, callback=callback)

        top_hits: TopHits
        for top_hits in hmmsearch_iter:
            for reported_hit in top_hits.reported:
                for domain in reported_hit.domains:
                    cds_idx = int(reported_hit.name.decode())
                    accession = domain.alignment.hmm_accession.decode()
                    score = domain.score
                    yield HSP(genes[cds_idx], accession, score)
