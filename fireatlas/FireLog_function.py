import logging
import os
from functools import wraps
from fireatlas.FireConsts import root_dir
from fireatlas import settings

_logger_configured = False

def get_logger(name):

    global _logger_configured

    logger = logging.getLogger(name)

    if not _logger_configured:
        logger.setLevel(logging.INFO)

        # create a console handler and set its level
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)

        # create a file handler as well
        fh = logging.FileHandler(os.path.join(root_dir, "running.log"))
        fh.setLevel(logging.INFO)

        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        ch.setFormatter(formatter)
        fh.setFormatter(formatter)

        logger.addHandler(ch)
        logger.addHandler(fh)

        # To avoid duplicate log messages when using `getLogger` with the same name,
        # prevent further propagation of messages to the root logger
        logger.propagate = False

        logger.info("logger initialized!")

        # Make file handler an attribute so you can change it later
        #logger.fh = fh
        
        _logger_configured = True

    return logger

def update_fh(logger, dirpath):
    logger.removeHandler(logger.fh)
    newpath = logger.default_file_path.replace(root_dir, dirpath)
    os.makedirs(dirpath, exist_ok=True)
    logger.fh = logging.FileHandler(newpath)
    logger.create_handler(logger.fh)
    return logger


def reset_fh(logger):
    dirpath = os.path.dirname(os.path.join(root_dir, "running.log"))
    return update_fh(logger, dirpath)

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
