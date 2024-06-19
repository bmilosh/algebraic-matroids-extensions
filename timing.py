from datetime import timedelta
from time import localtime, strftime, time


def secondsToStr(elapsed=None):
    if elapsed is None:
        return strftime("%Y-%m-%d %H:%M:%S", localtime())
    else:
        return str(timedelta(seconds=elapsed))

def log(s, elapsed=None):
    time_in_str = secondsToStr()
    line = "=" * (len(s) + len(time_in_str) + 3)
    print(line)
    print(secondsToStr(), '-', s)
    if elapsed:
        print("Elapsed time:", elapsed)
    print(line)
    print()

def endlog(start):
    end = time()
    elapsed = end-start
    log("End Program", secondsToStr(elapsed))
