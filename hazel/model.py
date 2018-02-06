from .atmosphere import Hazel_atmosphere, SIR_atmosphere, Parametric_atmosphere
from collections import OrderedDict
from . import pyhazel
from .spectrum import Spectrum
import numpy as np
import matplotlib.pyplot as pl
from ipdb import set_trace as stop
__all__ = ['Model']

class Model(object):
    def __init__(self):
        
        self.photospheres = []
        self.chromospheres = []
        self.chromospheres_order = []
        self.stray = []
        self.parametric = []
        self.spectrum = []

        # Initialize pyhazel
        pyhazel._init()

    def add_atmospheres(self, atmospheres, order=None):
        """
        Add a new photosphere to the model. They will be added with filling factors

        Parameters
        ----------
        photospheres : list 
            list of SIR_atmospheres that will be added with a filling factor

        chromospheres : list 
            list of lists of Hazel_atmospheres. Each element of the outer list refers to  that will be added with a filling factor
        
        Returns
        -------
        None

        """

        assert order in ['vertical', 'horizontal', None]

        n_atm = len(atmospheres)

        types = [atm.type for atm in atmospheres]

        atm_type = types[0]

        if (not all(x==atm_type for x in types)):
            raise ValueError("Added atmospheres are a mixture of photosphere and chromospheres and are not compatible.")
        else:

            if (atm_type == 'chromosphere'):
                self.chromospheres.append(atmospheres)
                if (len(atmospheres) > 1):
                    self.chromospheres_order.append(order)
                else:
                    self.chromospheres_order.append(None)

            if (atm_type == 'photosphere'):
                self.photospheres.append(atmospheres)
                

        # for atm in self.chromospheres
        
    # def add_chromosphere(self, order=0):
    #     """
    #     Add a new chromosphere to the model. The order is used to add chromospheres ones
    #     after the other

    #     Parameters
    #     ----------
    #     chromosphere : Hazel_atmosphere
    #         A Hazel atmosphere
    #     order : int (default: 0)
    #         Order, starting from 0
        
    #     Returns
    #     -------
    #     None

    #     """
    #     chromosphere = Hazel_atmosphere()
    #     self.chromospheres[order].append(chromosphere)

    def add_stray(self):
        """
        Add a new stray component

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        chromosphere = Hazel_atmosphere()
        self.chromospheres[order].append(chromosphere)

    def add_parametric(self):
        """
        Add a new parametric component

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        parametric = Parametric_atmosphere()
        self.parametric.append(parametric)

    def initialize_ff(self):
        """
        Normalize all filling factors so that they add to one to avoid later problems
        """
        
        ff = 0.0
        for atm in self.photospheres:
            ff += atm.ff

        for atm in self.photospheres:
            atm.ff /= ff


        ff = 0.0
        for atm in self.chromospheres:
            ff += atm.ff

        for atm in self.chromospheres:
            atm.ff /= ff

    def add_spectrum(self, spectra):
        """
        Add wavelength axis for the synthesis

        Parameters
        ----------
        spectra : float
            Wavelength axis vector
        
        Returns
        -------
        None
        """

        for sp in spectra:
            self.spectrum.append(sp)

    def list_wavelength_axis(self):
        """
        List all wavelength axis for the synthesis

        Parameters
        ----------
        None
        
        Returns
        -------
        None
        """

        n_axis = len(self.spectrum)
        print("N. wavelength axis : {0}".format(n_axis))
        for n, spectrum in enumerate(self.spectrum):
            print(" - Axis {0} : min={1} -> max={2}".format(n, np.min(spectrum.wavelength_axis), np.max(spectrum.wavelength_axis)))

    def synthesize(self):
        stokes = None
        if (self.photospheres is not None):
            for i, components in enumerate(self.photospheres):                
                for j, atm in enumerate(components):
                    stokes = atm.synthesize()

                    ind_low, ind_top = atm.wvl_range
                    atm.spectrum.stokes[:, ind_low:ind_top+1] = stokes[1:,:]

                    pl.plot(atm.spectrum.wavelength_axis, atm.spectrum.stokes[0,:])

        if (len(self.chromospheres) != 0):
            for i, components in enumerate(self.chromospheres):
                
                if (self.chromospheres_order[i] in ['vertical', None]):
                    for j, atm in enumerate(components):
                        if (j == 0):
                            stokes = atm.synthesize()
                        else:
                            stokes = atm.synthesize(stokes)

                        ind_low, ind_top = atm.wvl_range
                        atm.spectrum.stokes[:, ind_low:ind_top] = stokes

                        pl.plot(atm.spectrum.wavelength_axis, atm.spectrum.stokes[0,:])

        pl.show()
                                         
        
        # if (self.stray is not None):
            # for atm in self.stray:
                # atm.synthesize()

        # if (self.parametric is not None):
            # for atm in self.parametric:
                # atm.synthesize()