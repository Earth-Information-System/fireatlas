from functools import wraps
from time import time

def timed(f):
    from .FireLog import logger

    @wraps(f)
    def wrap(*args, **kwargs):    
        # capture the start time
        t_start = time()

        # run the function
        result = f(*args, **kwargs)

        # capture the end time
        t_end = time()

        # get the time difference
        t_diff = (t_end - t_start)

        # humanize the time difference
        if t_diff > 60:
            took = f"{t_diff / 60:.2f} min"
        elif t_diff > 1:
            took = f"{t_diff:.2f} sec"
        else:
            took = f"{t_diff * 1000:.2f} ms"
        
        # log the time that the function took
        logger.info(f"func:{f.__name__} took: {took}")
        return result
    return wrap