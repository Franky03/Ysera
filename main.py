import logging
import multiprocessing
from new_ysera import myfunction, mythread, lerDados
import os
import numpy as np
import time
from multiprocessing import Process,freeze_support, set_start_method,Pool,Manager
import sys
import multiprocessing as mp

if __name__ == '__main__':
    try:
        set_start_method('spawn')
    except RuntimeError:
        pass
    start_time = time.time()
    print(np.__version__)


    PROJECT_HOME = os.path.dirname(os.path.realpath(__file__))
    name = 'file_30.pdb'
    params = {}
    string = ""
    string2 = ""
    string3 = ""
    string4 = ""
    New, params = lerDados(name, params)

    print(f"New: {New}")
    print(f"Params: {params}")

    mythread(New, params, 0, 3515, "arquivo1", string)
    mythread(New, params, 3516, 7031, "arquivo2", string2)
    mythread(New, params, 7031, 10545, "arquivo3", string3)
    mythread(New, params, 10547, 14060, "arquivo4", string4)

    #Pool method not working

    # pool = Pool(processes=4)
    # r1 = pool.apply_async(mythread, [New, params, 0, 3515, "arquivo1", string])

    # r2 = pool.apply_async(mythread, [New, params, 3516, 7031, "arquivo2", string2])

    # r3 = pool.apply_async(mythread, [New, params, 7031, 10545, "arquivo3", string3])

    # r4 = pool.apply_async(mythread, [New, params, 10547, 14060, "arquivo4", string4])

    # pool.close()
    # pool.join()

    print("---%s seconds ---" % (time.time() - start_time))

#Linhas do arquivo 2 = 177, 614, 1289