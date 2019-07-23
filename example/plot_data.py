import numpy as np
import h5py

import glob

files = glob.glob('snap*')
Nfiles = len(files)

dat = [None] * Nfiles

for i in range(Nfiles):
    with h5py.File('snap_%03d.hdf5'%i,'r') as f:
        dat[i] = f.get('PartType1/Coordinates')[:]


dat = np.vstack(dat)

import matplotlib.pyplot as plt

mask = np.random.rand(int(len(dat))) < 1
plt.scatter(dat[mask,0], dat[mask,1], s=1)
plt.show()


### read grid
w = np.loadtxt('grid_weights.txt')
w = w.reshape((100,100,100))

w_collapse = w.sum(axis=2)

plt.imshow(w_collapse)
plt.show()


