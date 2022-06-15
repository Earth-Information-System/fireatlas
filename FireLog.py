""" FireLog
Module containing all logging info used in this project

Usage: from FireLog import logger

"""
from FireConsts import dirdata
import logging

# get logger
logger = logging.getLogger(__name__)

# set level
logger.setLevel(logging.INFO)

# define the console handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# define the file handler
fh = logging.FileHandler(dirdata+'running.log')  # the logfile is stored in the dirpjdata directory
fh.setLevel(logging.INFO)

# format the console and file handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
fh.setFormatter(formatter)

# add the handlers to logger
# logger.addHandler(ch)   # comment this to stop screen output
logger.addHandler(fh)   # comment this to stop log file recording
