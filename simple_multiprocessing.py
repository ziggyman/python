from multiprocessing import Pool
import numpy as np
import random # my pixels have lots of data in the center so I want to avoid getting one Thread having to deal with all the large files

def process(processNumber):
    print('processNumber = ',processNumber)

if __name__ == '__main__':
    p = Pool(processes=16)
    processNumbers = np.arange(100)
    random.shuffle(processNumbers)
    p.map(process, processNumbers)
    p.close()
