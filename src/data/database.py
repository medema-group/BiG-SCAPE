#!/usr/bin/env python
# vim: set fileencoding=utf-8 :
#
# Copyright (C) 2019 Satria A. Kautsar
# Wageningen University & Research
# Bioinformatics Group
"""bigslice.modules.data.database

Common classes and functions to work with the SQLite3 database
"""

from os import path
from shutil import move
import re
import sqlite3
from time import time


class Database:
    """Wrapper class to do manipulation on an SQLite3 database file"""

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()

    def __init__(self, db_path: str, use_memory: bool=False,
                 for_query_mode: bool=False):
        """db_path: path to sqlite3 database file
        CAUTION: this should not be subjected to
        multiple processes (at least for the insert
        function"""

        self._db_path = db_path
        self._insert_queues = []
        self._insert_queues_index = {}
        self._last_indexes = {}
        self._connection = None
        self._use_memory = use_memory

        if for_query_mode:
            # truncated database object for query mode
            sql_schema = open(path.join(path.dirname(
                path.abspath(__file__)), "schema_query_mode.sql"), "r").read()
            if path.exists(self._db_path):
                print(("database file {} exists! this is not "
                       "what it was made for.").format(self._db_path))
            else:
                self._connection = sqlite3.connect(self._db_path)
                db_cur = self._connection.cursor()
                db_cur.executescript(sql_schema)
            for row in self.select("sqlite_sequence", "WHERE 1"):
                self._last_indexes[row["name"]] = row["seq"]
            return

        # get schema information
        sql_schema = open(path.join(path.dirname(
            path.abspath(__file__)), "schema.sql"), "r").read()
        self.schema_ver = re.search(
            r"\n-- schema ver\.: (?P<ver>\d+\.\d+\.\d+)", sql_schema).group("ver")

        if path.exists(self._db_path):
            if (self._use_memory):
                with sqlite3.connect(self._db_path) as db_con:
                    query = "".join(line for line in db_con.iterdump())
                    self._connection = sqlite3.connect(":memory:")
                    self._connection.executescript(query)
            else:
                self._connection = sqlite3.connect(self._db_path)
            # check if existing one have the same schema version
            db_schema_ver = self.select("schema", "WHERE 1")[0]["ver"]
            if db_schema_ver != self.schema_ver:
                raise Exception(
                    "SQLite3 database exists but contains different schema " +
                    "version ({} rather than {}), exiting!".format(
                        db_schema_ver, self.schema_ver))
        else:
            # create new database
            if (self._use_memory):
                self._connection = sqlite3.connect(":memory:")
            else:
                self._connection = sqlite3.connect(self._db_path)
            db_cur = self._connection.cursor()
            db_cur.executescript(sql_schema)

            # fetch last_indexes from db
            for row in self.select("sqlite_sequence", "WHERE 1"):
                self._last_indexes[row["name"]] = row["seq"]

            # load bs_class rows into a dictionary for quick searching
            bs_classes_id = {}
            for row in self.select(
                "chem_class",
                "WHERE 1"
            ):
                bs_classes_id[row["name"]] = row["id"]
                if row["name"] == "Unknown":
                    # insert chem_subclass: Unknown-unknown
                    self.insert("chem_subclass",
                                {
                                    "name": "unknown",
                                    "class_id": row["id"]
                                })
                    self.commit_inserts()

            # load chem_class_map.tsv
            with open(path.join(path.dirname(path.abspath(__file__)),
                                "chem_class_map.tsv"), "r") as tsv:
                tsv.readline()  # skip header
                for line in tsv:
                    src_class, src, bs_class, bs_subclass = \
                        line.rstrip().split("\t")
                    if bs_subclass == "":  # TODO: enforce strict
                        bs_subclass = "unknown"
                    if bs_class == "":  # TODO: enforce strict
                        bs_class = "Unknown"
                    bs_class_id = bs_classes_id[bs_class]

                    existing = self.select(
                        "chem_subclass",
                        "WHERE name=? AND class_id=?",
                        parameters=(bs_subclass, bs_class_id)
                    )
                    if existing:
                        assert len(existing) == 1
                        bs_subclass_id = existing[0]["id"]
                    else:
                        bs_subclass_id = self.insert(
                            "chem_subclass",
                            {
                                "name": bs_subclass,
                                "class_id": bs_class_id
                            }
                        )

                    self.insert("chem_subclass_map", {
                        "class_source": src_class,
                        "type_source": src,
                        "subclass_id": bs_subclass_id
                    })
                self.commit_inserts()

        # fetch last_indexes from db
        for row in self.select("sqlite_sequence", "WHERE 1"):
            self._last_indexes[row["name"]] = row["seq"]

    def close(self):
        if self._use_memory:
            self.dump_db_file()
        self._connection.close()

    def dump_db_file(self):
        if self._use_memory:
            start = time()
            print("Dumping in-memory database content into " +
                  self._db_path + "...", end=" ", flush=True)
            if path.exists(self._db_path):
                move(self._db_path, self._db_path + ".bak")
            with sqlite3.connect(self._db_path) as out_db:
                query = "".join([line for line in self._connection.iterdump()])
                out_db.executescript(query)
            print("{0:.4f}s".format(time() - start))
        else:
            raise(Exception("not an in-memory database"))

    def select(self, table: str, clause: str,
               parameters: tuple = None, props: list = [],
               distinct: bool = False,
               as_tuples: bool = False):
        """execute a SELECT ... FROM ... WHERE"""

        if len(props) < 1:
            props_string = "*"
        else:
            props_string = ",".join(props)

        if distinct:
            props_string = "DISTINCT " + props_string

        sql = "SELECT {} FROM {} {}".format(
            props_string,
            table,
            clause
        )

        def dict_factory(cursor, row):
            """see https://docs.python.org/2/library/
            sqlite3.html#sqlite3.Connection.row_factory"""
            d = {}
            for idx, col in enumerate(cursor.description):
                d[col[0]] = row[idx]
            return d

        orig_row_factory = self._connection.row_factory
        if not as_tuples:
            self._connection.row_factory = dict_factory
        db_cur = self._connection.cursor()

        if parameters:
            results = db_cur.execute(sql, parameters).fetchall()
        else:
            results = db_cur.execute(sql).fetchall()

        self._connection.row_factory = orig_row_factory
        return results

    def update(self, table: str, data: dict, clause: str,
               parameters: tuple = ()):
        """execute an UPDATE ... SET ... WHERE ..."""

        set_clause = ""
        set_params = []
        for k, v in data.items():
            if set_clause:
                set_clause += ","
            set_clause += k + "=?"
            set_params.append(v)
        set_params.extend(parameters)

        sql = "UPDATE {} SET {} {}".format(
            table,
            set_clause,
            clause
        )

        db_cur = self._connection.cursor()
        db_cur.execute(sql, tuple(set_params))
        self._connection.commit()

        return db_cur.rowcount

    def insert(self, table: str, data: dict):
        """execute an INSERT INTO ... VALUES ...
        !!don't use the returned IDs unless you are sure
        that the INSERTs all will be successful and there
        is no other sources tinkering with the database"""

        new_id = self._last_indexes.get(table, 0) + 1
        self._last_indexes[table] = new_id

        self._insert_queues.append((table, data, new_id))
        if table not in self._insert_queues_index:
            self._insert_queues_index[table] = []
        self._insert_queues_index[table].append(len(self._insert_queues) - 1)
        return new_id

    def get_pending_id(self, table: str, query: dict):
        """try to look for a match in the insert buffer,
        and returns the pending id
        !!don't use this unless you are sure
        that the INSERTs all will be successful and there
        is no other sources tinkering with the database"""

        ids = []
        for buffer_idx in self._insert_queues_index.get(table, []):
            _, buffer_data, pending_id = self._insert_queues[buffer_idx]
            match = True
            for key in query:
                if key not in buffer_data or \
                        buffer_data[key] != query[key]:
                    match = False
                    break
            if match:
                ids.append(pending_id)
        return ids

    def commit_inserts(self):
        """perform actual commit for the insert queue"""

        db_cur = self._connection.cursor()
        for table, data, _ in self._insert_queues:
            keys = []
            values = []
            for key, value in data.items():
                if key == "id":  # can't have this around!
                    raise Exception("Don't specify id for INSERTs!")
                keys.append(key)
                values.append(value)
            sql = "INSERT INTO {}({}) VALUES ({})".format(
                table,
                ",".join(keys),
                ",".join(["?" for i in range(len(values))])
            )
            db_cur.execute(sql, tuple(values))
        self._connection.commit()
        self._insert_queues = []
        self._insert_queues_index = {}

        # sanity check
        for row in self.select("sqlite_sequence", "WHERE 1"):
            if self._last_indexes[row["name"]] != row["seq"]:
                raise Exception("buffered indexes no longer " +
                                "in sync with sqlite_sequence, database " +
                                "might be corrupted!!!")
