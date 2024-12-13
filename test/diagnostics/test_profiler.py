"""module containing a test for the profiler"""

# from python
from unittest import TestCase
import psutil
import os
import time

# from this module
from big_scape.diagnostics.profiler import calc_cpu_percent


class TestProfiler(TestCase):
    """Test class for the profiler functions"""

    def test_cpu_percent_within_range(self):
        """Tests if the calc_cpu_percent method outputs a value between 0-100"""
        process = psutil.Process(os.getpid())
        start_time = process.cpu_times().user
        update_interval = 1
        time.sleep(update_interval)

        cpu_percent = calc_cpu_percent(process, start_time, update_interval)

        self.assertTrue(0 <= cpu_percent <= 100)
