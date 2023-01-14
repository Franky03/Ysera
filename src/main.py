from ysera_pandas import *
from thread import Thread
from multiprocessing import Pool
import time
from multiprocessing import Manager
import glob
import os


def run_ysera(file):
    print('Loading Data...')
    af = AromaticsFormat(filename=file)
    aromatic_array, aromatic_normals, invalids, \
        total, total_dist = af.get_data()
    print(total_dist)
    print('Data Loaded')
    tr = Thread(aromatic_array, aromatic_normals, invalids, total, total_dist)
    size = len(total)
    print('Calculating Bonds...')
    if size > 3000:
        p = int(size/4)
        r = size % 4
        d = 0
        parts = []
        for i in range(4):
            parts.append(p + d)
            d += p
        parts[3] = parts[3] + r
        manager = Manager()
        output_dict = manager.dict()
        pool = Pool(processes=4)
        pool.apply_async(tr.run, (0, parts[0], 'ysera_1', output_dict))
        pool.apply_async(tr.run, (parts[0], parts[1], 'ysera_2', output_dict))
        pool.apply_async(tr.run, (parts[1], parts[2], 'ysera_3', output_dict))
        pool.apply_async(tr.run, (parts[2], parts[3], 'ysera_4', output_dict))
        pool.close()
        pool.join()
        output_name = filename.replace('.pdb', '')
        final_text = output_dict.values()
        final_text = ''.join(final_text)
        with open(f'output/{output_name}.txt', 'w') as f:
            f.write(final_text)


filename = 'file_21.pdb'
if __name__ == '__main__':
    start = time.time()
    run_ysera(filename)
    print(f"---{(time.time() - start)} seconds ---")
