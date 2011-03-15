import time

def print_timing(func):
    def wrapper(*arg,**kwargs):
        t1 = time.time()
        res = func(*arg,**kwargs)
        t2 = time.time()
        print '%s took %0.5g s' % (func.func_name, (t2-t1))
        return res
    return wrapper

