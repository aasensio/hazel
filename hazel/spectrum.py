import numpy as np

__all__ = ['Spectrum']

class Spectrum(object):
    def __init__(self, wvl=None, weights=None, observed_file=None):
        
        self.wavelength_axis = None
        self.stokes = None

        if (wvl is not None):
            self.add_spectrum(wvl)

        if (weights is not None):
            self.add_weights(weights)

        if (observed_file is not None):
            self.add_observed_file(observed_file)

    def add_spectrum(self, wvl):
        self.wavelength_axis = wvl
        self.stokes = np.zeros((4,len(wvl)))

    def add_weights(self, weights):
        self.wavelength_weights = weights

    def add_observed_file(self, observed_file):
        self.observed_file = observed_file