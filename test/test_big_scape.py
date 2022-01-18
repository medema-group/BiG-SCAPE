import argparse
import os
import sys

from unittest import TestCase
import unittest

from src import big_scape
from src import utility

class test_options:
    label: str

class TestRunBase(TestCase):
    test_run: big_scape.Run

    def setUp(self):
        self.test_run = big_scape.Run()
        self.test_run.run_mode = "Test"

        # test options object
        self.test_run.options = test_options()
        self.test_run.options.label = "Test"

    def test_start(self):
        self.test_run.start(True)

        self.assertIsNotNone(self.test_run.start_time)
    
    def test_end(self):
        self.test_run.start(True)
        self.test_run.end()

        self.assertIsNotNone(self.test_run.run_data["end_time"])
