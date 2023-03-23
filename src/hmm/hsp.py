"""Contains class to describe domain hits and alignments"""


# from other modules
from src.genbank import CDS


class HSP:
    """Describes a CDS - Domain relationship"""

    def __init__(self, cds: CDS, domain: str, score) -> None:
        self.cds = cds
        self.domain = domain
        self.score = score

    def save(self):
        """Saves this object to a database"""
        pass

    def __repr__(self) -> str:
        return ",".join(
            map(
                str,
                [
                    self.cds.nt_start,
                    self.cds.nt_stop,
                    self.cds.strand,
                    self.domain,
                    self.score,
                ],
            )
        )

    def __eq__(self, __o: object) -> bool:
        if not isinstance(__o, HSP):
            return False

        return all(
            [
                __o.cds.nt_start == self.cds.nt_start,
                __o.cds.nt_stop == self.cds.nt_stop,
                __o.domain == self.domain,
                __o.score == self.score,
            ]
        )


class HSPAlignment:
    """Describes a CDS - Domain alignment"""

    def __init__(self) -> None:
        pass

    def save(self):
        """Saves this object to a database"""
        pass
