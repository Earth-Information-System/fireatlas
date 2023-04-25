import functools
import subprocess
from FireLog import logger


def free_profiler(target_func):
    @functools.wraps(target_func)
    def wrapper(*args, **kwargs):
        try:
            result = subprocess.check_output(['bash', '-c', 'free -t -m -h'])
            logger.info(result)
        except:
            pass
        return target_func(*args, **kwargs)
    return wrapper