"""
This module provides basic loggers
"""
import logging
import os

def basic_logger(log2stdout=logging.INFO, log2file=None, name='el_default',
                 logFileName="log.txt", purge_old=True):
    _logger = logging.getLogger(name)
    _logger.setLevel(log2stdout)
    msg_str = '%(asctime)s::%(name)s::%(lineno)s::%(levelname)s - %(message)s'
    msg_str = '%(funcName)s::%(filename)s:%(lineno)s::%(levelname)s - %(message)s'
    date_format = "%H:%M:%S"
    formatter = logging.Formatter(msg_str, datefmt=date_format)
    if log2file is not None:
        if purge_old and os.path.exists(logFileName):
            os.unlink(logFileName)
        fh = logging.FileHandler(logFileName)
        fh.setLevel(log2file)
        fh.setFormatter(formatter)
        _logger.addHandler(fh)
    if log2stdout is not None:
        ch = logging.StreamHandler()
        ch.setLevel(log2stdout)
        ch.setFormatter(formatter)
        _logger.addHandler(ch)
    return _logger

