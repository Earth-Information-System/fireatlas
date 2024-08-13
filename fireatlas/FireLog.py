import logging
import os
from fireatlas import settings
from fireatlas.FireConsts import root_dir

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

    def create_handler(self, handler, level):
        handler.setLevel(level)
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        self.addHandler(handler)

    def set_fh(self, path):
        self.fh = logging.FileHandler(path)
        self.create_handler(self.fh, logging.INFO)

    def configure_logger(self):
        self.setLevel(logging.INFO)

        # create a console handler and set its level
        ch = logging.StreamHandler()
        self.create_handler(ch, logging.INFO)

        # create a file handler as well
        self.set_fh(self.default_file_path)

        # To avoid duplicate log messages when using `getLogger` with the same name,
        # prevent further propagation of messages to the root logger
        self.propagate = False # not sure if this does anything actually
        self.info("logger initialized!")

    def update_fh_dir(self, dirpath):
        self.removeHandler(self.fh)
        newpath = self.default_file_path.replace(root_dir, dirpath)
        os.makedirs(dirpath, exist_ok=True)
        self.set_fh(newpath)
    
    def reset_fh(self):
        self.update_fh_dir(os.path.dirname(self.default_file_path)) 


# Register the custom logger class with the logging module
logging.setLoggerClass(CustomLogger)

# create logger instance
logger = logging.getLogger(__name__)      