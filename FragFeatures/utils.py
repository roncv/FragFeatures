"""
Various utility functions for FragFeatures.
"""
import time


# Decorator to time functions
def timeit(func):
    """
    Decorator to time functions.
    """
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        print(f"function {func.__name__} took {end - start:.2f}s to run.")
        return result
    return wrapper
