import numpy as np
import matplotlib.pyplot as pl
import hazel
import h5py
from ipdb import set_trace as stop

# Single invert
mod = hazel.Model('conf_single_invert.ini')
mod.invert()