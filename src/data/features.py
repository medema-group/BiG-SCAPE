#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2019 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group
"""bigslice.modules.data.features

Handle manipulation and storage of features
for now, also do feature extraction

todo: separate features extraction method
todo: use parallel processing
todo: use Dask
"""

from os import path
from .database import Database
from multiprocessing import Pool
import numpy as np
import pandas as pd


class Features:
    """Represents a features entry in the database"""

    def __init__(self, properties: dict):
        self.bgc_id = properties["bgc_id"]
        self.hmm_id = properties["hmm_id"]
        self.value = properties["value"]

    def save(self, database: Database):
        """commits bgc_features extraction data"""
        existing = database.select(
            "bgc_features",
            "WHERE bgc_id=? and hmm_id=?",
            parameters=(self.bgc_id, self.hmm_id),
            props=["bgc_id"]
        )
        if existing:
            # for now, this should not get called
            raise Exception("not_implemented")
        else:
            # save
            database.insert(
                "bgc_features",
                {
                    "bgc_id": self.bgc_id,
                    "hmm_id": self.hmm_id,
                    "value": self.value
                }
            )

    @staticmethod
    def extract(bgc_ids: int,
                database: Database):
        """ Extract features from domain hits
        """

        # prepare features extraction
        hmm_ids, parent_hmm_ids, hmm_names = map(
            tuple, list(zip(*database.select(
                "hmm LEFT JOIN subpfam ON hmm.id=subpfam.hmm_id",
                "",
                props=["hmm.id", "subpfam.parent_hmm_id", "hmm.name"],
                as_tuples=True
            ))))

        # for referencing the rows
        hmm_idx = {value: idx for idx, value in enumerate(hmm_ids)}
        subpfam_ids = {}
        for i, hmm_id in enumerate(hmm_ids):
            parent_hmm_id = parent_hmm_ids[i]
            if parent_hmm_id:
                if parent_hmm_id not in subpfam_ids:
                    subpfam_ids[parent_hmm_id] = []
                subpfam_ids[parent_hmm_id].append(hmm_id)

        # fetch features
        features = []
        hsps = {bgc_id: {} for bgc_id in bgc_ids}
        for bgc_id, cds_id, hmm_id, bitscore in database.select(
            "hsp,cds",
            "WHERE hsp.cds_id = cds.id" +
            " AND cds.bgc_id in (" + ",".join(map(str, bgc_ids)) + ")" +
            " AND hsp.hmm_id in (" + ",".join(map(str, hmm_ids)) + ")",
            props=["cds.bgc_id", "cds.id", "hsp.hmm_id",
                   "CAST(bitscore as INTEGER)"],
            as_tuples=True
        ):
            if cds_id not in hsps[bgc_id]:
                hsps[bgc_id][cds_id] = {}
            if hmm_id not in hsps[bgc_id][cds_id]:
                hsps[bgc_id][cds_id][hmm_id] = []
            hsps[bgc_id][cds_id][hmm_id].append(bitscore)

        for bgc_id in hsps:
            biosyn_present = set()  # per-bgc
            subpfam_bitscores = {}  # per-hmm, per cds-region
            for cds_id in hsps[bgc_id]:
                for hmm_id in hsps[bgc_id][cds_id]:
                    bitscore = max(hsps[bgc_id][cds_id][hmm_id])
                    parent_hmm_id = parent_hmm_ids[hmm_idx[hmm_id]]
                    if not parent_hmm_id:  # biosyn
                        biosyn_present.add(hmm_id)
                    else:  # subpfam
                        subpfam_bitscores[hmm_id] = int(
                            np.max(hsps[bgc_id][cds_id][hmm_id]))

            for hmm_id in biosyn_present:
                features.append(Features({
                    "bgc_id": bgc_id,
                    "hmm_id": hmm_id,
                    "value": 255
                }))

            for hmm_id, bitscore in subpfam_bitscores.items():
                features.append(Features({
                    "bgc_id": bgc_id,
                    "hmm_id": hmm_id,
                    "value": bitscore
                }))

        return features
