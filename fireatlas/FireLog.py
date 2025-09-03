import logging
import os
from fireatlas import settings
import sys 
import platform 
import subprocess
import datetime as dt 

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
        fh = logging.FileHandler(settings.LOG_FILEPATH)
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


def write_run_metadata():
    """
    Writes a text file to settings.ENV_META_FILEPATH that captures the output of 
    pip freeze and basic platform information. 
    """
    t = dt.datetime.utcnow().isoformat() + "Z"
    python_version = sys.version.replace("\n", " ")
    platform_info = platform.platform() 

    meta = [
        f"# Environment metadata generated: {t}",
        f"# Python: {python_version}", 
        f"# Platform: {platform_info}", 
        "#" * 60
    ]

    freeze_output = subprocess.run(
        ["conda", "env", "export"], 
        capture_output=True, 
        text=True, 
        check=True
    )

    content = "\n".join(meta) + "\n" + freeze_output.stdout

    with open(settings.ENV_META_FILEPATH, "w", encoding="utf-8") as f:
        f.write(content) 
    
    return settings.ENV_META_FILEPATH 
