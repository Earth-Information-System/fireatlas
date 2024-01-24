from functools import wraps
from time import time
from FireLog import logger

def timed(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        t_start = time()
        result = f(*args, **kwargs)
        t_end = time()
        t_diff = (t_end - t_start)
        if t_diff > 60:
            took = f"{t_diff / 60:.2f} min"
        elif t_diff > 1:
            took = f"{t_diff:.2f} sec"
        else:
            took = f"{t_diff * 1000:.2f} ms"
        logger.info(f"func:{f.__name__} took: {took}")
        return result
    return wrap