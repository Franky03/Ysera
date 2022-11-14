from ysera_pandas import AromaticsFormat
from thread import Thread
from multiprocessing import Pool
import time


def run_ysera(file):
    print('Loading Data...')
    af = AromaticsFormat(filename=file)
    aromatic_array, aromatic_normals, invalids, \
        total, total_dist = af.get_data()
    print('Data Loaded')
    p = int(len(total)/4)
    r = len(total) % 4
    d = 0
    parts = []
    for i in range(4):
        parts.append(p + d)
        d += p
    parts[3] = parts[3] + r
    tr = Thread(aromatic_array, aromatic_normals, invalids, total, total_dist)
    print('Calculating Bonds...')
    pool = Pool(processes=4)
    pool.apply_async(tr.run, (0, parts[0], 'ysera_1'))
    pool.apply_async(tr.run, (parts[0], parts[1], 'ysera_2'))
    pool.apply_async(tr.run, (parts[1], parts[2], 'ysera_3'))
    pool.apply_async(tr.run, (parts[2], parts[3], 'ysera_4'))
    pool.close()
    pool.join()


filename = 'file_30.pdb'
if __name__ == '__main__':
    start = time.time()
    run_ysera(filename)
    print(f"---{(time.time() - start)} seconds ---")
# 358s arquivo completo
