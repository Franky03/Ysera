from new_ysera import mythread, lerDados
import os
import numpy as np
import time
import asyncio

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
        loop= asyncio.get_event_loop()

        r1= mythread(New, params, 0, 3515, f"arquivo1_{name[5:7]}", string)
        r2= mythread(New, params, 3516, 7031, f"arquivo2_{name[5:7]}", string2)
        r3= mythread(New, params, 7031, 10545, f"arquivo3_{name[5:7]}", string3)
        r4= mythread(New, params, 10547, 14060, f"arquivo4_{name[5:7]}", string4)

        all_rs= asyncio.gather(r1,r2,r3,r4)

        Result= loop.run_until_complete(all_rs)

        loop.close()
        with open("output2.0/Ysera2.0.txt", "wb") as outfile:
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
        path = os.path.dirname(os.path.realpath(__file__))
        os.remove(path + '/output2.0/ysera_1.txt')
        os.remove(path + '/output2.0/ysera_2.txt')
        os.remove(path + '/output2.0/ysera_3.txt')
        os.remove(path + '/output2.0/ysera_4.txt')
        print(Result)
    except Exception as e:
        print(e)

    print("---%s seconds ---" % (time.time() - start_time))

#Linhas do arquivo 2 = 177, 614, 1289