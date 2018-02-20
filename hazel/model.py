from hazel.atmosphere import Hazel_atmosphere, SIR_atmosphere, Parametric_atmosphere
from hazel.configuration import Configuration
from collections import OrderedDict
from hazel.codes import hazel_code
from hazel.spectrum import Spectrum
import hazel.util
import numpy as np
import matplotlib.pyplot as pl
from pathlib import Path
import h5py
from ipdb import set_trace as stop
__all__ = ['Model']

class Model(object):
    def __init__(self, config=None, verbose=None):
        
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

            if (verbose is not None):
                self.verbose = verbose
            else:
                self.verbose = bool(self.configuration.config_dict['working mode']['verbose'])

            self.use_configuration(self.configuration.config_dict)

        # Initialize pyhazel
        hazel_code._init()

    def use_configuration(self, config_dict):

        # Deal with the spectral regions        
        tmp = config_dict['spectral regions']

        self.output_file = config_dict['working mode']['output file']

        spec = {}

        for key, value in tmp.items():
            if (self.verbose):
                print('Adding spectral region {0}'.format(value['name']))

            if (value['wavelength file'] == 'None'):
                if ('lower, upper, n. wavelengths' in value):
                    axis = value['lower, upper, n. wavelengths']
                    wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
            else:
                if (self.verbose):
                    print('  - Reading wavelength axis from {0}'.format(value['wavelength file']))
                wvl = np.loadtxt(value['wavelength file'])

            if (value['wavelength weight file'] == 'None'):
                if (self.verbose):
                    print('  - Setting all weights to 1')
                weights = np.ones(len(wvl))
            else:
                if (self.verbose):
                    print('  - Reading wavelength weights from {0}'.format(value['wavelength weight file']))
                weights = np.loadtxt(value['wavelength weight file'])

            if (value['observations file'] == 'None'):
                if (self.verbose):
                    print('  - Not using observations')
                obs_file = None
            else:
                if (self.verbose):
                    print('  - Using observations from {0}'.format(value['observations file']))
                obs_file = value['observations file']
                        
            spec[value['name']] = Spectrum(wvl, weights, obs_file)
        
        self.spectrum = spec
        
        # Deal with the atmospheres
        tmp = config_dict['atmospheres']

        self.atmospheres = {}

        append = False
        index_chromosphere = 0
        index_photosphere = 0

        if (self.verbose):
            print('Adding atmospheres')

        for key, value in tmp.items():
            if ('photosphere' in key):
                if (self.verbose):
                    print('  - New available photosphere : {0}'.format(value['name']))
                
                self.atmospheres[value['name']] = SIR_atmosphere()
                lines = [int(k) for k in list(value['spectral lines'])]
                wvl_range = [float(k) for k in value['wavelength range']]
                            
                self.atmospheres[value['name']].add_active_line(index=index_photosphere, lines=lines, spectrum=spec[value['spectral region']], 
                    wvl_range=np.array(wvl_range), append=append)

                my_file = Path(value['reference atmospheric model'])
                if (not my_file.exists()):
                    raise FileExistsError("Input file for atmosphere {0} does not exist.".format(value['name']))

                self.atmospheres[value['name']].load_reference_model(value['reference atmospheric model'], self.verbose)

                if (self.atmospheres[value['name']].model_type == '3d'):                    
                    f = h5py.File(self.atmospheres[value['name']].model_file, 'r')
                    self.atmospheres[value['name']].n_pixel, _, _ = f['model'].shape
                    f.close()
                    
                for k, v in value['nodes'].items():
                    for k2, v2 in self.atmospheres[value['name']].parameters.items():
                        if (k.lower() == k2.lower()):
                            tmp, cycles = hazel.util._extract_parameter_cycles(v)
                            self.atmospheres[value['name']].parameters[k2] = tmp
                            self.atmospheres[value['name']].cycles[k2] = tmp
                
                index_photosphere += 1
            
            if ('chromosphere' in key):
                if (self.verbose):
                    print('  - New available chromosphere : {0}'.format(value['name']))
                
                self.atmospheres[value['name']] = Hazel_atmosphere()
                
                wvl_range = [float(k) for k in value['wavelength range']]

                self.atmospheres[value['name']].add_active_line(index=index_chromosphere, line=value['line'], spectrum=spec[value['spectral region']], 
                    wvl_range=np.array(wvl_range))

                my_file = Path(value['reference atmospheric model'])
                if (not my_file.exists()):
                    raise FileExistsError("Input file for atmosphere {0} does not exist.".format(value['name']))

                self.atmospheres[value['name']].load_reference_model(value['reference atmospheric model'], self.verbose)

                if (self.atmospheres[value['name']].model_type == '3d'):
                    f = h5py.File(self.atmospheres[value['name']].model_file, 'r')
                    self.atmospheres[value['name']].n_pixel, _ = f['model'].shape
                    f.close()
                
                # Set values of parameters
                self.atmospheres[value['name']].parameters['h'] = float(value['height'])
                for k, v in value['parameters'].items():
                    for k2, v2 in self.atmospheres[value['name']].parameters.items():
                        if (k.lower() == k2.lower()):                            
                            self.atmospheres[value['name']].cycles[k2] = v

                index_chromosphere += 1
            

        tmp = config_dict['topology']
        for key, value in tmp.items():
            self.add_atmospheres(value)
            
        to_remove = []        
        for k, v in self.atmospheres.items():
            if (not v.active):
                to_remove.append(k)
                if (self.verbose):
                    print('  - Atmosphere {0} is not used. Deleted.'.format(k))
                
        for k in to_remove:
            self.atmospheres.pop(k)

        self.initialize_ff()

        # Check that number of pixels is the same for all read files
        n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
        all_equal = all(x == n_pixels[0] for x in n_pixels)
        if (not all_equal):
            for k, v in self.atmospheres.items():
                print('{0} -> {1}'.format(k, v.n_pixel))
            raise Exception("Files with model atmospheres do not contain the same number of pixels")
        else:
            if (self.verbose):
                print('  - Number of pixels to read : {0}'.format(n_pixels[0]))
            self.n_pixels = n_pixels[0]
                    
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
            print('  - Adding topology : {0}'.format(atmosphere_order))

        vertical_order = atmosphere_order.split('->')        
        order = []
        for k in vertical_order:
            name = k.strip().replace('(','').replace(')','').split('+')
            
            tmp = []
            for n in name:
                tmp.append(self.atmospheres[n])
                self.atmospheres[n].active = True

            order.append(tmp)            

        self.order_atmospheres.append(order)
                
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
        if (self.verbose):
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
                                                            
                    if (n > 0):
                        ind_low, ind_top = atm.wvl_range
                        stokes_out = atm.spectrum.stokes[:, ind_low:ind_top]
                                                                            
                    if (k == 0):
                        stokes = atm.ff * atm.synthesize(stokes_out)
                    else:                        
                        stokes += atm.ff * atm.synthesize(stokes_out)
                                
                ind_low, ind_top = atm.wvl_range
                
                atm.spectrum.stokes[:, ind_low:ind_top+1] = stokes

    def plot_stokes(self):        
        for atmospheres in self.order_atmospheres:
            f, ax = pl.subplots(nrows=2, ncols=2)
            ax = ax.flatten()
            for i in range(4):
                ax[i].plot(atmospheres[0][-1].spectrum.stokes[i,:])
