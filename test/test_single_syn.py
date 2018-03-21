import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop


# Single

mod = hazel.Model('conf_single.ini')
mod.synthesize()
f, ax = pl.subplots(nrows=2, ncols=2)
ax = ax.flatten()
for i in range(4):
    ax[i].plot(mod.spectrum['spec1'].stokes[i,:])
pl.show()
pl.pause(0.001)

np.savetxt('10830_stokes.1d', mod.spectrum['spec1'].stokes.T + 1e-4 * np.random.randn(150,4), header='lambda SI SQ SU SV')
stop()
