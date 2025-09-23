# Python script to store python outputs to file and terminal
# Author: S.W.Hoh, University of York, 2024

import sys
from typing import Literal

Mode = Literal['w', 'w+', 'a']


class log2file():
    """
    Logs to console and file using the python logging module

    Usage:
    >>> log = log2file("output.log", "w", shell=True)
    >>> sys.stdout = sys.stderr = log
    >>> print("stuff")
    """
    def __init__(self, filename: str, mode: Mode = 'w', shell: bool = True):
        """Set up logger

        Args:
            filename (str): Log file path
            mode (Mode, optional): Write mode, 'w', 'w+', 'a'. Defaults to 'w'.
            shell (bool, optional): Output to shell/terminal. Defaults to True.
        """
        #self.level = level
        #self.linebuf = ""
        #self.stream = stream
        self.to_term = shell
        self.file = open(filename, mode)
        self.stdout = sys.__stdout__
        self.stderr = sys.__stderr__

    def write(self, message):
        """Writes output to file and terminal if log2file(..., shell=True)

        Args:
            message (str): Message to be written
        """
        if self.to_term:
            self.stdout.write(message)
        self.file.write(message)

    def flush(self):
        """Flush std output or error
        """
        if self.to_term:
            self.stdout.flush()
            self.stderr.flush()
        self.file.flush()
        
    def close(self):
        """Close file
        """
        self.file.close()
    
    @staticmethod
    def reset_stdouterr_to_sys():
        sys.stdout = sys.__stdout__
        sys.stderr = sys.__stderr__

'''def set_logger(
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
    logger = logging.getLogger(filename)
    handler = logging.FileHandler(filename)
    logger.addHandler(handler)
    stdout_level *= 10
    stderr_level *= 10
    sys.stdout = log2file(logger, stdout_level, orig_stdout) 
    sys.stderr = log2file(logger, stderr_level, orig_stderr) 
    
    #log = logging.getLogger(__name__)
    #stdout_level *= 10
    #stderr_level *= 10
    #sys.stdout = log2file(log, stdout_level, orig_stdout)
    #sys.stderr = log2file(log, stderr_level, orig_stderr)
'''