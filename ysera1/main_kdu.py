from ysera_pandas import AromaticsFormat
from ysera1.thread import Thread
from multiprocessing import Pool
import time
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

        pool = Pool(processes=4)
        pool.apply_async(tr.run, (0, parts[0], 'ysera_1'))
        pool.apply_async(tr.run, (parts[0], parts[1], 'ysera_2'))
        pool.apply_async(tr.run, (parts[1], parts[2], 'ysera_3'))
        pool.apply_async(tr.run, (parts[2], parts[3], 'ysera_4'))
        pool.close()
        pool.join()
        read_files = glob.glob("output_2.0/*.txt")
        with open("output_2.0/Ysera2.0.txt", "wb") as outfile:
            for f in read_files:
                with open(f, "rb") as infile:
                    outfile.write(infile.read())
        path = os.path.dirname(os.path.realpath(__file__))
        os.remove(path + '/output_2.0/ysera_1.txt')
        try:
            os.remove(path + '/output_2.0/ysera_2.txt')
        except IOError:
            pass
        try:
            os.remove(path + '/output_2.0/ysera_3.txt')
        except IOError:
            pass
        try:
            os.remove(path + '/output_2.0/ysera_4.txt')
        except IOError:
            pass
    else:
        tr.run(0, len(total), 'Ysera2.0')


filename = '3og7.pdb'
if __name__ == '__main__':
    start = time.time()
    run_ysera(filename)
    print(f"---{(time.time() - start)} seconds ---")
