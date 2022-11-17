from ysera_pandas import AromaticsFormat
from thread import Thread
from multiprocessing import Pool, Manager
import time
import glob
import os


def run_ysera(file):
    print('Loading Data...')
    af = AromaticsFormat(filename=file)
    aromatic_array, aromatic_normals, invalids, \
        total, total_dist = af.get_data()
    print('Data Loaded')
    manager = Manager()
    exclusions = manager.list()
    dist_old = manager.list()
    tr = Thread(aromatic_array, aromatic_normals, invalids, total, total_dist, exclusions, dist_old)
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
        pool = Pool(processes=4)
        pool.apply_async(tr.run, (0, parts[0], 'ysera_1'))
        time.sleep(5)
        pool.apply_async(tr.run, (parts[0], parts[1], 'ysera_2'))
        time.sleep(5)
        pool.apply_async(tr.run, (parts[1], parts[2], 'ysera_3'))
        time.sleep(5)
        pool.apply_async(tr.run, (parts[2], parts[3], 'ysera_4'))
        pool.close()
        pool.join()
        read_files = glob.glob("output2.0/*.txt")
        with open("output2.0/Ysera2.0.txt", "wb") as outfile:
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
        path = os.path.dirname(os.path.realpath(__file__))
        # os.remove(path + '/output2.0/ysera_1.txt')
        # os.remove(path + '/output2.0/ysera_2.txt')
        # os.remove(path + '/output2.0/ysera_3.txt')
        # os.remove(path + '/output2.0/ysera_4.txt')
    else:
        tr.run(0, len(total), 'Ysera2.0')


filename = 'file_30.pdb'
if __name__ == '__main__':
    start = time.time()
    run_ysera(filename)
    print(f"---{(time.time() - start)} seconds ---")
