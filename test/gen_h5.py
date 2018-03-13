import numpy as np
import h5py

n_pixel = 50

model = np.loadtxt('model_photosphere.1d', dtype=np.float32, skiprows=1)

nz, ncol = model.shape
model_3d = np.zeros((n_pixel,nz,ncol), dtype=np.float32)

for i in range(n_pixel):
    model[:,1] += 100.0
    model_3d[i,:,:] = model
    
f = h5py.File('model_photosphere.h5', 'w')
db = f.create_dataset('model', model_3d.shape, dtype=np.float32)
db[:] = model_3d
f.close()

model = np.loadtxt('model_chromosphere.1d', skiprows=1)

model_3d = np.zeros((n_pixel,8), dtype=np.float64)
for i in range(n_pixel):
    model_3d[i,:] = model

f = h5py.File('model_chromosphere.h5', 'w')
db = f.create_dataset('model', model_3d.shape, dtype=np.float64)
db[:] = model_3d
f.close()

model = np.loadtxt('model_chromosphere2.1d', skiprows=1)

model_3d = np.zeros((n_pixel,8), dtype=np.float64)
for i in range(n_pixel):
    model[5] += 2.0
    model_3d[i,:] = model
    
f = h5py.File('model_chromosphere2.h5', 'w')
db = f.create_dataset('model', model_3d.shape, dtype=np.float64)
db[:] = model_3d
f.close()