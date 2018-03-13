import numpy as np
from hazel.io import Generic_observed_file

__all__ = ['Spectrum']

class Spectrum(object):
    def __init__(self, wvl=None, weights=None, observed_file=None, name=None):
        
        self.wavelength_axis = None
        self.stokes = None
        self.stokes_perturbed = None

        if (wvl is not None):
            self.add_spectrum(wvl)

        if (weights is not None):
            self.add_weights(weights)

        if (observed_file is not None):
            self.add_observed_file(observed_file)

        if (name is not None):
            self.add_name(name)

    def add_spectrum(self, wvl):
        self.wavelength_axis = wvl
        self.stokes = np.zeros((4,len(wvl)))
        self.stokes_perturbed = np.zeros((4,len(wvl)))

    def add_weights(self, weights):
        self.wavelength_weights = weights

    def add_observed_file(self, observed_file):        
        self.observed_handle = Generic_observed_file(observed_file)

    def add_name(self, name):
        self.name = name