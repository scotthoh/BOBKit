# Python script to store python outputs to file
# Author: S.W.Hoh, University of York, 2024

import sys
import logging


class log2file(object):
    """
    Logs to console and file using the python logging module

    Usage:
    >>> logging.basicConfig(level=logging.DEBUG,
                            filename="testout.log",
                            filemode="a",
                            format="%(message)s")
    >>> log = logging.getLogger(__name__)
    >>> sys.stdout = log2file(log, logging.INFO, orig_stdout)
    >>> sys.stderr = log2file(log, logging.INFO, orig_stdout)
    >>> print("stuff")

    """

    def __init__(self, logger, level, stream):
        """
        Set up things
        """
        self.logger = logger
        self.level = level
        self.linebuf = ""
        self.stream = stream

    def write(self, data):
        """
        writes output
        """
        for line in data.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

        self.stream.write(data)

    def flush(self):
        pass


def set_logger(
    filename: str = "output_log.txt",
    filemode: str = "a",
    stdout_level: int = 1,
    stderr_level: int = 4,
):
    """Set up logger to writes stdout and stderr from python to file,
    stdout/stderr from C++ are not logged into this file.

    Args:
        filename (str, optional): Name of output file to write to.
        filemode (str, optional): w for write, a for append.
                                  Defaults to a for append.
        stdout_level (int, optional): Numeric value for stdout logging level,
                                      refer logging level in python.
                                      Defaults to 1.
        stderr_level (int, optional): Numeric value for stderr logging level,
                                      refer logging level in python.
                                      Defaults to 4.

    """
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    # auto set to append the file if filemode not set
    if filemode not in ["a", "w"]:
        filemode = "a"
    logging.basicConfig(
        level=logging.DEBUG,
        filename=filename,
        filemode=filemode,
        format="%(message)s",
    )
    log = logging.getLogger(__name__)
    stdout_level *= 10
    stderr_level *= 10
    sys.stdout = log2file(log, stdout_level, orig_stdout)
    sys.stderr = log2file(log, stderr_level, orig_stderr)
