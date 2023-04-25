import functools
import subprocess
from FireLog import logger


def free_profiler(target_func):
    """a quick and dirty `free -m` memory profiler
    b/c backgrounding bash commands from the main process
    wasn't working. Note this decorator will slow down your code

    :param target_func:
    :return:
    """
    @functools.wraps(target_func)
    def wrapper(*args, **kwargs):
        try:
            results = subprocess.check_output(['bash', '-c', 'free -t -m -h'])
            split_lines = results.split(b'\n')
            if not split_lines: return
            for result_line in split_lines:
                logger.info(result_line)
        except Exception as e:
            logger.exception(e)
        return target_func(*args, **kwargs)
    return wrapper
