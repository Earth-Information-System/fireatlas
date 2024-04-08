import logging
import os

_logger_configured = False

def get_logger(name):
    from fireatlas.FireConsts import root_dir

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
        _logger_configured = True

    return logger

logger = get_logger(__name__)
