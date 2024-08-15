import logging
import os
from fireatlas import settings
from fireatlas.FireConsts import root_dir
from functools import wraps

# singleton ensures logger is initiated only once
class Singleton(type):
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
    
class CustomLogger(logging.Logger, metaclass=Singleton):
    def __init__(self, name):
        super().__init__(name)
        self.default_file_path = os.path.join(root_dir, settings.LOG_FILENAME)
        self.configure_logger()

    def create_handler(self, handler):
        handler.setLevel(logging.INFO)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.addHandler(handler)

    def configure_logger(self):
        self.setLevel(logging.INFO)

        # Remove any existing NullHandler instances
        self.handlers = [h for h in self.handlers if not isinstance(h, logging.NullHandler)]


        # create a console handler and set its level
        ch = logging.StreamHandler()
        self.create_handler(ch)

        # create a file handler as well
        self.fh = logging.FileHandler(self.default_file_path)
        self.create_handler(self.fh)

        # To avoid duplicate log messages when using `getLogger` with the same name,
        # prevent further propagation of messages to the root logger
        self.propagate = False # not sure if this does anything actually
        self.info("logger initialized!")

    def update_fh_dir(self, dirpath):
        self.removeHandler(self.fh)
        newpath = self.default_file_path.replace(root_dir, dirpath)
        os.makedirs(dirpath, exist_ok=True)
        self.fh = logging.FileHandler(newpath)
        self.create_handler(self.fh)
    
    def reset_fh(self):
        self.update_fh_dir(os.path.dirname(self.default_file_path)) 

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
                logger.update_fh_dir(all_dir_func(tst, region, location = None))

            result = f(*args, **kwargs)
            
                # reset logging path to default
            if settings.LOG_SUBDIR and tst and region:
                logger.reset_fh()
                
            return result
        return wrapper
    return decorator

# Register the custom logger class with the logging module
logging.setLoggerClass(CustomLogger)

# create logger instance
logger = logging.getLogger(__name__)      