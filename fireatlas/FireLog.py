import logging
import os
from fireatlas import settings
from fireatlas.FireConsts import root_dir
from functools import wraps 

DEFAULT_FILE_PATH = os.path.join(root_dir, settings.LOG_FILENAME)

_logger_configured = False

def create_handler(logger, handler):
    handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger

def get_logger(name):

    global _logger_configured

    logger = logging.getLogger(name)

    if not _logger_configured:
        logger.setLevel(logging.INFO)

        # create a console handler and set its level
        ch = logging.StreamHandler()
        create_handler(logger, ch)

        # create a file handler as well
        fh = logging.FileHandler(DEFAULT_FILE_PATH)
        create_handler(logger, fh)
        
        # add file handler as an attribute in order to possibly update
        logger.fh = fh

        # To avoid duplicate log messages when using `getLogger` with the same name,
        # prevent further propagation of messages to the root logger
        logger.propagate = False
        logger.info("logger initialized!")

        _logger_configured = True

    return logger


def update_fh(logger, dirpath):
    # remove the old one
    logger.removeHandler(logger.fh)
    newpath = DEFAULT_FILE_PATH.replace(os.path.dirname(DEFAULT_FILE_PATH), dirpath)
    os.makedirs(dirpath, exist_ok=True)

    # create new file handler
    logger.fh = logging.FileHandler(newpath)
    create_handler(logger, logger.fh)

    return logger


def reset_fh(logger):
    return update_fh(logger, os.path.dirname(DEFAULT_FILE_PATH))

# a decorator to apply to certain fuctions in order to change the file handler path
def logger_subdir(all_dir_func, tst_pos, region_pos):
    def decorator(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            # Extract tst and region from kwargs or args
            tst = kwargs.get('tst', args[tst_pos] if len(args) > tst_pos else None)
            region = kwargs.get('region', args[region_pos] if len(args) > region_pos else None)

            # configure logger to put logs into subdirectory
            if settings.LOG_SUBDIR and tst and region:
                update_fh(logger, all_dir_func(tst, region, location = None))

            result = f(*args, **kwargs)
            
                # reset logging path to default
            if settings.LOG_SUBDIR and tst and region:
                reset_fh(logger)
                
            return result
        return wrapper
    return decorator

logger = get_logger(__name__)
