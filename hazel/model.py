from hazel.atmosphere import Hazel_atmosphere, SIR_atmosphere, Parametric_atmosphere
from hazel.configuration import Configuration
from collections import OrderedDict
from hazel.codes import hazel_code
from hazel.spectrum import Spectrum
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
        hazel_code._init()

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

        append = False
        index_chromosphere = 0
        index_photosphere = 0

        for key, value in tmp.items():
            if ('photosphere' in key):
                if (self.verbose):
                    print("Adding photosphere : {0}".format(value['name']))
                self.atmospheres[value['name']] = SIR_atmosphere()
                lines = [int(k) for k in list(value['spectral lines'])]
                wvl_range = [float(k) for k in value['wavelength range']]
                            
                self.atmospheres[value['name']].add_active_line(index=index_photosphere, lines=lines, spectrum=spec[value['spectral region']], 
                    wvl_range=np.array(wvl_range), append=append)
                self.atmospheres[value['name']].load_model(value['input model'])

                index_photosphere += 1
            
            if ('chromosphere' in key):
                if (self.verbose):
                    print("Adding chromosphere : {0}".format(value['name']))
                    self.atmospheres[value['name']] = Hazel_atmosphere()
                    
                    wvl_range = [float(k) for k in value['wavelength range']]

                    self.atmospheres[value['name']].add_active_line(index=index_chromosphere, line=value['line'], spectrum=spec[value['spectral region']], wvl_range=np.array(wvl_range))

                    # Set values of parameters
                    self.atmospheres[value['name']].parameters['h'] = float(value['height'])
                    for k, v in value['parameters'].items():
                        for k2, v2 in self.atmospheres[value['name']].parameters.items():
                            if (k.lower() == k2.lower()):
                                self.atmospheres[value['name']].parameters[k2] = float(v)

                    index_chromosphere += 1

        tmp = config_dict['topology']
        for key, value in tmp.items():
            self.add_atmospheres(value)


        self.initialize_ff()
        # self.adapt_wavelength_photosphere()

    def adapt_wavelength_photosphere(self):
        n_lambda = []
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                n_lambda.append(len(v.wvl_axis))

        shift = 0
        ind = 0
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                v.wvl_range += shift
                shift += n_lambda[ind]
                ind += 1
                        

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
        if (self.verbose):
            print("Adding atmosphere : {0}".format(atmosphere_order))

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
        """
        Synthesize all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        stokes = None
        stokes_out = None

        for atmospheres in self.order_atmospheres:
            
            for n, order in enumerate(atmospheres):
                                                                
                for k, atm in enumerate(order):

                    print(k, atm)
                                                            
                    if (n > 0):
                        ind_low, ind_top = atm.wvl_range
                        stokes_out = atm.spectrum.stokes[:, ind_low:ind_top]
                                                                            
                    if (k == 0):
                        stokes = atm.ff * atm.synthesize(stokes_out)
                    else:                        
                        stokes += atm.ff * atm.synthesize(stokes_out)
                                
                ind_low, ind_top = atm.wvl_range
                
                atm.spectrum.stokes[:, ind_low:ind_top+1] = stokes