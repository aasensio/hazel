from .atmosphere import Hazel_atmosphere, SIR_atmosphere, Parametric_atmosphere
from .configuration import Configuration
from collections import OrderedDict
from . import pyhazel
from .spectrum import Spectrum
import numpy as np
import matplotlib.pyplot as pl
from ipdb import set_trace as stop
__all__ = ['Model']

class Model(object):
    def __init__(self, config=None):
        
        self.photospheres = []
        self.chromospheres = []
        self.chromospheres_order = []
        self.atmospheres = {}
        self.order_atmospheres = []
        self.stray = []
        self.parametric = []
        self.spectrum = []
        self.configuration = None

        if (config is not None):
            self.configuration = Configuration(config)
            self.use_configuration(self.configuration.config_dict)

        # Initialize pyhazel
        pyhazel._init()

    def use_configuration(self, config_dict):

        # Deal with the spectral regions        
        tmp = config_dict['spectral regions']

        self.verbose = bool(config_dict['working mode']['verbose'])

        spec = {}

        for key, value in tmp.items():
            axis = value['lower, upper, n. wavelengths']
            wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
            spec[value['name']] = Spectrum(wvl)
        
        spec_list = [v for k, v in spec.items()]

        self.add_spectrum(spec)

        # Deal with the atmospheres
        tmp = config_dict['atmospheres']

        self.atmospheres = {}

        for key, value in tmp.items():
            if ('photosphere' in key):
                if (self.verbose):
                    print("Adding photosphere : {0}".format(value['name']))
                self.atmospheres[value['name']] = SIR_atmosphere()
                lines = [int(k) for k in list(value['spectral lines'])]
                wvl_range = [float(k) for k in value['wavelength range']]
                            
                self.atmospheres[value['name']].add_active_line(lines=lines, spectrum=spec[value['spectral region']], wvl_range=np.array(wvl_range))
                self.atmospheres[value['name']].load_model(value['input model'])
                
            
            if ('chromosphere' in key):
                if (self.verbose):
                    print("Adding chromosphere : {0}".format(value['name']))
                    self.atmospheres[value['name']] = Hazel_atmosphere()
                    
                    wvl_range = [float(k) for k in value['wavelength range']]

                    self.atmospheres[value['name']].add_active_line(line=value['line'], spectrum=spec[value['spectral region']], wvl_range=np.array(wvl_range))

                    # Set values of parameters
                    self.atmospheres[value['name']].parameters['h'] = float(value['height'])
                    for k, v in value['parameters'].items():
                        for k2, v2 in self.atmospheres[value['name']].parameters.items():
                            if (k.lower() == k2.lower()):
                                self.atmospheres[value['name']].parameters[k2] = float(v)

        tmp = config_dict['topology']
        for key, value in tmp.items():
            self.add_atmospheres(value)


        self.initialize_ff()


    def add_atmospheres(self, atmosphere_order):
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

        # Transform the order to a list of lists        
        vertical_order = atmosphere_order.split('->')        
        order = []
        for k in vertical_order:
            name = k.strip().replace('(','').replace(')','').split('+')
            
            tmp = []
            for n in name:
                tmp.append(self.atmospheres[n])

            order.append(tmp)

        self.order_atmospheres.append(order)        

        # assert order in ['vertical', 'horizontal', None]

        # n_atm = len(atmospheres)

        # types = [atm.type for atm in atmospheres]

        # atm_type = types[0]

        # if (not all(x==atm_type for x in types)):
        #     raise ValueError("Added atmospheres are a mixture of photosphere and chromospheres and are not compatible.")
        # else:

        #     if (atm_type == 'chromosphere'):
        #         self.chromospheres.append(atmospheres)
        #         if (len(atmospheres) > 1):
        #             self.chromospheres_order.append(order)
        #         else:
        #             self.chromospheres_order.append(None)

        #     if (atm_type == 'photosphere'):
        #         self.photospheres.append(atmospheres)
                
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

        for atmospheres in self.order_atmospheres:
            for order in atmospheres:
                total_ff = 0.0
                for atm in order:
                    total_ff += atm.ff

                for atm in order:                    
                    atm.ff /= total_ff                


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
        stokes_out = None

        for atmospheres in self.order_atmospheres:
            
            for n, order in enumerate(atmospheres):
                                                                
                for k, atm in enumerate(order):
                                                            
                    if (n > 0):
                        ind_low, ind_top = atm.wvl_range
                        stokes_out = atm.spectrum.stokes[:, ind_low:ind_top]
                                                                            
                    if (k == 0):
                        stokes = atm.ff * atm.synthesize(stokes_out)
                    else:                        
                        stokes += atm.ff * atm.synthesize(stokes_out)
                                
                ind_low, ind_top = atm.wvl_range
                
                atm.spectrum.stokes[:, ind_low:ind_top+1] = stokes

                pl.plot(atm.spectrum.wavelength_axis, atm.spectrum.stokes[0,:])
        stop()
                
                
                    

                    

        # # for component in self.order_atmospheres:
            
        # if (self.photospheres is not None):
        #     for i, components in enumerate(self.photospheres):                
        #         for j, atm in enumerate(components):
        #             stokes = atm.synthesize()

        #             ind_low, ind_top = atm.wvl_range
        #             atm.spectrum.stokes[:, ind_low:ind_top+1] = stokes[1:,:]

        #             pl.plot(atm.spectrum.wavelength_axis, atm.spectrum.stokes[0,:])

        # if (len(self.chromospheres) != 0):
        #     for i, components in enumerate(self.chromospheres):
                
        #         if (self.chromospheres_order[i] in ['vertical', None]):
        #             for j, atm in enumerate(components):
        #                 if (j == 0):
        #                     stokes = atm.synthesize()
        #                 else:
        #                     stokes = atm.synthesize(stokes)

        #                 ind_low, ind_top = atm.wvl_range
        #                 atm.spectrum.stokes[:, ind_low:ind_top] = stokes

        #                 pl.plot(atm.spectrum.wavelength_axis, atm.spectrum.stokes[0,:])

        pl.show()
                                         
        
        # if (self.stray is not None):
            # for atm in self.stray:
                # atm.synthesize()

        # if (self.parametric is not None):
            # for atm in self.parametric:
                # atm.synthesize()