import numpy as np

__all__ = ['Spectrum']

class Spectrum(object):
    def __init__(self, wvl=None):
        
        self.wavelength_axis = None
        self.stokes = None

        if (wvl is not None):
            self.add_spectrum(wvl)

    def add_spectrum(self, wvl):
        
        self.wavelength_axis = wvl
        self.stokes = np.zeros((4,len(wvl)))
