import numpy as np
import hazel.h5py_pickle as h5py
from astropy.io import fits
import os
from ipdb import set_trace as stop

__all__ = ['Generic_observed_file', 'Generic_hazel_file', 'Generic_SIR_file', 'Generic_parametric_file', 'Generic_stray_file']

class Generic_observed_file(object):

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            tmp = np.loadtxt(self.filename, skiprows=1).T
            return tmp[0:4,:], tmp[4:,:]

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler.close()
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_stray_file(object):

    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            tmp = np.loadtxt(self.filename, skiprows=1)
            return tmp

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler.close()
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_parametric_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            return np.loadtxt(self.filename, skiprows=1)

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler.close()
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_hazel_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return

        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            return np.loadtxt(self.filename, skiprows=1)

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

        # if (self.extension == 'fits'):
        #     return self.handler[0]fits.open(self.filename, memmap=True)
        #     return

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler.close()
        

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _ = self.handler['model'].shape
            self.close()
            return tmp

class Generic_SIR_file(object):
    
    def __init__(self, filename):
        self.extension = os.path.splitext(filename)[1][1:]
        self.filename = filename

    def open(self):
        if (self.extension == '1d'):
            return
        if (self.extension == 'h5'):
            self.handler = h5py.File(self.filename, 'r')
            return
        if (self.extension == 'fits'):
            self.handler = fits.open(self.filename, memmap=True)
            return

    def read(self, pixel=None):
        if (self.extension == '1d'):
            f = open(self.filename, 'r')
            f.readline()
            ff = float(f.readline())
            f.close()
            return np.loadtxt(self.filename, skiprows=4), ff

        if (self.extension == 'h5'):
            return self.handler['model'][pixel,...]

    def close(self):
        if (self.extension == '1d'):
            return

        if (self.extension == 'h5'):
            self.handler.close()

    def get_npixel(self):
        if (self.extension == '1d'):
            return 1

        if (self.extension == 'h5'):
            self.open()
            tmp, _, _ = self.handler['model'].shape
            self.close()
            return tmp