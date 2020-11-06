import functools
import time
import numpy as np
import pandas as pd

from sth_simulation.helsim_RUN import STH_Simulation

def timer(func):
    """Print the runtime of the decorated function"""
    @functools.wraps(func)
    def wrapper_timer(*args, **kwargs):
        print(f"-> Running {func.__name__!r}")
        start_time = time.perf_counter()    # 1
        value = func(*args, **kwargs)
        end_time = time.perf_counter()      # 2
        run_time = end_time - start_time    # 3
        print(f"=> Finished {func.__name__!r} in {run_time:.4f} secs\n\n")
        return value
    return wrapper_timer

@timer
def test_1a():
    print("Run the simulations for 7 years.")
    STH_Simulation(paramFileName='AscarisParameters_moderate.txt',
                   demogName='KenyaKDHS',
                   MDAFilePath='files/Input_MDA_23Oct20.csv',
                   prevKKSACFilePath='files/OutputPrevKKSAC_STH_test_1a.csv',
                   prevMHISACFilePath='files/OutputPrevMHISAC_STH_test_1a.csv',
                   RkFilePath='files/InputRk_STH.csv',
                   nYears=7,
                   outputFrequency=12,
                   numReps=None,
                   SaveOutput=False,
                   OutSimFilePath=None,
                   InSimFilePath=None)
    print("results in files/OutputPrevKKSAC_STH_test_1a.csv and files/OutputPrevMHISAC_STH_test_1a.csv")

@timer
def test_1b():
    print("Rerun the same simulations also for 7 years.")
    STH_Simulation(paramFileName='AscarisParameters_moderate.txt',
                   demogName='KenyaKDHS',
                   MDAFilePath='files/Input_MDA_23Oct20.csv',
                   prevKKSACFilePath='files/OutputPrevKKSAC_STH_test_1b.csv',
                   prevMHISACFilePath='files/OutputPrevMHISAC_STH_test_1b.csv',
                   RkFilePath='files/InputRk_STH.csv',
                   nYears=7,
                   outputFrequency=12,
                   numReps=None,
                   SaveOutput=False,
                   OutSimFilePath=None,
                   InSimFilePath=None)
    print("results in files/OutputPrevKKSAC_STH_test_1b.csv and files/OutputPrevMHISAC_STH_test_1b.csv")

@timer
def test_2a():
    print("Rerun the same simulations for the first 2 years.")
    STH_Simulation(paramFileName='AscarisParameters_moderate.txt',
                   demogName='KenyaKDHS',
                   MDAFilePath='files/Input_MDA_23Oct20_part1.csv',
                   prevKKSACFilePath='files/OutputPrevKKSAC_STH_test_2a.csv',
                   prevMHISACFilePath='files/OutputPrevMHISAC_STH_test_2a.csv',
                   RkFilePath='files/InputRk_STH.csv',
                   nYears=2,
                   outputFrequency=12,
                   numReps=None,
                   SaveOutput=True,
                   OutSimFilePath='files/Output_STH_test_2a.p',
                   InSimFilePath=None)
    print("results in files/OutputPrevKKSAC_STH_test_2a.csv and files/OutputPrevMHISAC_STH_test_2a.csv")

@timer
def test_2b():
    print("Continue the same simulations for the remaining 5 years.")
    STH_Simulation(paramFileName='AscarisParameters_moderate.txt',
                   demogName='KenyaKDHS',
                   MDAFilePath='files/Input_MDA_23Oct20_part2.csv',
                   prevKKSACFilePath='files/OutputPrevKKSAC_STH_test_2b.csv',
                   prevMHISACFilePath='files/OutputPrevMHISAC_STH_test_2b.csv',
                   RkFilePath='files/InputRk_STH.csv',
                   nYears=5,
                   outputFrequency=12,
                   numReps=None,
                   SaveOutput=False,
                   OutSimFilePath=None,
                   InSimFilePath='files/Output_STH_test_2a.p')
    print("results in files/OutputPrevKKSAC_STH_test_2b.csv and files/OutputPrevMHISAC_STH_test_2b.csv")