import numpy as np
import h5py

Nfiles = 100
Npart = int(1e4)
boxl = 50.


for i in range(Nfiles):
    print('Processing file ',i)
    with h5py.File('snap_%03d.hdf5'%i,'w') as f:
        f.create_group('PartType1')

        ## Random uniform positions
        coods = np.random.rand(Npart,3) * boxl
        print('min, max coords = ',np.min(coods),np.max(coods))

        ## Sinusoidal velocities
        vels = np.cos(4.*np.pi*coods/boxl)
        

        f.create_dataset('PartType1/Coordinates',data=coods)
        f.create_dataset('PartType1/Velocity',data=vels)

        

