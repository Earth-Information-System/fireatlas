import logging
import os
from fireatlas import settings
from functools import wraps 

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
        fh = logging.FileHandler(settings.LOG_PATH)
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
    newpath = settings.LOG_PATH.replace(os.path.dirname(settings.LOG_PATH), dirpath)
    os.makedirs(dirpath, exist_ok=True)

    # create new file handler
    logger.fh = logging.FileHandler(newpath)
    create_handler(logger, logger.fh)

    return logger


def reset_fh(logger):
    return update_fh(logger, os.path.dirname(settings.LOG_PATH))

def logger_subdir(all_dir_func, tst_pos: int, region_pos: int):
    ''' a decorator to be applied to certain fuctions in order to change the file handler path.
    logging file will be located in the FEDSoutput subdirectory for a given region/year.
    Subdirectory path resets after function runs.
    
    Parameters
    ----------

    all_dir_func : `preprocess.all_dir`. Unfortunately this function needs to be called outside of 
    `FireLog` due to circular dependencies. `logger_subdir` is not designed to work with any other 
    function. 
    tst_pos: index position of the tst argument in whatever function this decorator is applied to
    region_pos: index position of the region argument in whatever function this decorator is applied to

    '''
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
