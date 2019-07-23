import numpy as np
import h5py

Nfiles = 100
Npart = int(1e3)
boxl = 100


for i in range(Nfiles):
    with h5py.File('snap_%03d.hdf5'%i,'w') as f:
        f.create_group('PartType1')

        ## Random noise
        coods = np.random.rand(Npart,3) * boxl

        ## Gaussian
        # coods = np.random.normal(loc=boxl/2,scale=boxl/4,size=(Npart,3))
        import scipy.stats as stats
        mu, cov = boxl/2,boxl/4
        a, b = (0 - mu) / cov, (boxl - mu) / cov

        X = stats.truncnorm(a,b,boxl/2,boxl/4)
        coods =X.rvs((Npart,3))

        f.create_dataset('PartType1/Coordinates',data=coods)

        

