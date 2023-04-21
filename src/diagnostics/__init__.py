"""Contains modules for diagnostics"""

from src.diagnostics.logger import init_logger, init_logger_file
from src.diagnostics.profiler import Profiler, calc_cpu_percent

__all__ = [
    "init_logger",
    "init_logger_file",
    "Profiler",
    "calc_cpu_percent",
    "get_stats",
    "child_start_cpu_time",
    "make_plots",
    "collect_consumption",
]
