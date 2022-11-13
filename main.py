from new_ysera import myfunction, mythread, lerDados
import os
import numpy as np
import time
import asyncio

async def main():

    L= await asyncio.gather(
        mythread(New, params, 0, 3515, f"arquivo1_{name[5:7]}", string),
        mythread(New, params, 3516, 7031, f"arquivo2_{name[5:7]}", string2),
        mythread(New, params, 7031, 10545, f"arquivo3_{name[5:7]}", string3),
        mythread(New, params, 10547, 14060, f"arquivo4_{name[5:7]}", string4)
    )

    print(L)
        

if __name__ == '__main__':

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

    try:
        asyncio.run(main())

    except Exception as e:
        print(e)

    print("---%s seconds ---" % (time.time() - start_time))

#Linhas do arquivo 2 = 177, 614, 1289