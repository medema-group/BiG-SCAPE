"""Contains modules for diagnostics"""

from src.diagnostics.logger import init_logger, init_logger_file
from src.diagnostics.profiler import Profiler

__all__ = ["init_logger", "init_logger_file", "Profiler"]
