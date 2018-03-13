from hazel.chromosphere import Hazel_atmosphere
from hazel.photosphere import SIR_atmosphere
from hazel.parametric import Parametric_atmosphere
from hazel.configuration import Configuration
from collections import OrderedDict
from hazel.codes import hazel_code, sir_code
from hazel.spectrum import Spectrum
from hazel.transforms import transformed_to_physical, physical_to_transformed
import hazel.util
import numpy as np
import matplotlib.pyplot as pl
import copy
from pathlib import Path

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

        self.epsilon = 1e-1

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

        self.working_mode = config_dict['working mode']['action']

        spec = {}

        for key, value in tmp.items():
            if (self.verbose):
                print('Adding spectral region {0}'.format(value['name']))

            if (value['wavelength file'] == 'None'):
                if ('lower, upper, n. wavelengths' in value):
                    axis = value['lower, upper, n. wavelengths']
                    wvl = np.linspace(float(axis[0]), float(axis[1]), int(axis[2]))
                    print('  - Using wavelength axis from {0} to {1} with {2} steps'.format(float(axis[0]), float(axis[1]), int(axis[2])))
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
                        
            spec[value['name']] = Spectrum(wvl=wvl, weights=weights, observed_file=obs_file, name=value['name'])
        
        self.spectrum = spec
        
        # Deal with the atmospheres
        tmp = config_dict['atmospheres']

        self.atmospheres = {}

        if (self.verbose):
            print('Adding atmospheres')

        for key, value in tmp.items():
            
            if ('photosphere' in key):
                if (self.verbose):
                    print('  - New available photosphere : {0}'.format(value['name']))

                self.add_photosphere(value)
                                                            
            if ('chromosphere' in key):
                if (self.verbose):
                    print('  - New available chromosphere : {0}'.format(value['name']))

                self.add_chromosphere(value)
                                            
            if ('parametric' in key):
                if (self.verbose):
                    print('  - New available parametric : {0}'.format(value['name']))

                self.add_parametric(value)
                            

        tmp = config_dict['topology']
        if (self.verbose):
            print("Adding topologies") 
        for key, value in tmp.items():
            self.add_topology(value)

        # Remove unused atmospheres defined in the configuration file and not in the topology
        if (self.verbose):
            print("Removing unused atmospheres")
        self.remove_unused_atmosphere()
        

        # Calculate indices for atmospheres
        index_chromosphere = 1
        index_photosphere = 1
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                v.index = index_photosphere
                index_photosphere += 1
            if (v.type == 'chromosphere'):
                v.index = index_chromosphere
                index_chromosphere += 1

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
                print('Number of pixels to read : {0}'.format(n_pixels[0]))
            self.n_pixels = n_pixels[0]
        

        # Check that the number of pixels from all observations (in case of inversion) is the same
        # Check also if they are equal to those of the models
        # n_pixels = [v.n_pixel for k, v in self.atmospheres.items()]
        # all_equal = all(x == n_pixels[0] for x in n_pixels)

        # Check that the number of cycles is the same for all atmospheres (in case of inversion)
        if (self.working_mode == 'inversion'):
            cycles = []
            for k, v in self.atmospheres.items():
                for k2, v2 in v.cycles.items():                    
                    cycles.append(len(v2))

            all_equal = all(x == cycles[0] for x in cycles)
            if (not all_equal):
                raise Exception("Number of cycles in the nodes of active atmospheres is not always the same")
            else:
                self.n_cycles = cycles[0]

        self.init_sir()

    def add_photosphere(self, atmosphere):
        """
        Programmatically add a photosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength range', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = SIR_atmosphere()
        lines = [int(k) for k in list(atm['spectral lines'])]
        wvl_range = [float(k) for k in atm['wavelength range']]
                    
        self.atmospheres[atm['name']].add_active_line(lines=lines, spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        my_file = Path(atm['reference atmospheric model'])
        if (not my_file.exists()):
            raise FileExistsError("Input file for atmosphere {0} does not exist.".format(atm['name']))

        self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

        if (self.atmospheres[atm['name']].model_type == '3d'):
            self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
            
        for k, v in atm['nodes'].items():
            for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                if (k.lower() == k2.lower()):
                    self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)


    def add_chromosphere(self, atmosphere):
        """
        Programmatically add a chromosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength range', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)
        
        self.atmospheres[atm['name']] = Hazel_atmosphere()
                
        wvl_range = [float(k) for k in atm['wavelength range']]

        self.atmospheres[atm['name']].add_active_line(line=atm['line'], spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        if (self.verbose):
            print("    * Adding line : {0}".format(atm['line']))

        my_file = Path(atm['reference atmospheric model'])
        if (not my_file.exists()):
            raise FileExistsError("Input file for atmosphere {0} does not exist.".format(atm['name']))

        self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

        if (self.atmospheres[atm['name']].model_type == '3d'):
            self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        self.atmospheres[atm['name']].height = float(atm['height'])
        for k, v in atm['nodes'].items():
            for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                if (k.lower() == k2.lower()):                            
                    self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

    def add_parametric(self, atmosphere):
        """
        Programmatically add a parametric atmosphere

        Parameters
        ----------
        atmosphere : dict
            Dictionary containing the following data
            'Name', 'Spectral region', 'Height', 'Line', 'Wavelength range', 'Reference atmospheric model',
            'Ranges', 'Nodes'

        Returns
        -------
        None
        """

        # Make sure that all keys of the input dictionary are in lower case
        # This is irrelevant if a configuration file is used because this has been
        # already done
        atm = hazel.util.lower_dict_keys(atmosphere)

        self.atmospheres[atm['name']] = Parametric_atmosphere()
                
        wvl_range = [float(k) for k in atm['wavelength range']]

        self.atmospheres[atm['name']].add_active_line(spectrum=self.spectrum[atm['spectral region']], 
            wvl_range=np.array(wvl_range))

        my_file = Path(atm['reference atmospheric model'])
        if (not my_file.exists()):
            raise FileExistsError("Input file for atmosphere {0} does not exist.".format(atm['name']))

        self.atmospheres[atm['name']].load_reference_model(atm['reference atmospheric model'], self.verbose)

        if (self.atmospheres[atm['name']].model_type == '3d'):
            self.atmospheres[atm['name']].n_pixel = self.atmospheres[atm['name']].model_handler.get_npixel()
        
        # Set values of parameters
        for k, v in atm['nodes'].items():
            for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                if (k.lower() == k2.lower()):                            
                    self.atmospheres[atm['name']].cycles[k2] = hazel.util.toint(v)

        if ('ranges' in atm):
            for k, v in atm['ranges'].items():
                for k2, v2 in self.atmospheres[atm['name']].parameters.items():
                    if (k.lower() == k2.lower()):
                        if (v == 'None'):
                            self.atmospheres[atm['name']].ranges[k2] = None
                        else:
                            self.atmospheres[atm['name']].ranges[k2] = hazel.util.tofloat(v)

    def remove_unused_atmosphere(self):
        """
        Remove unused atmospheres
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        
        to_remove = []        
        for k, v in self.atmospheres.items():
            if (not v.active):
                to_remove.append(k)
                if (self.verbose):
                    print('  - Atmosphere {0} deleted.'.format(k))
                
        for k in to_remove:
            self.atmospheres.pop(k)
                    
    def init_sir(self):
        """
        Initialize SIR for this synthesis
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
    
        """
        for k, v in self.atmospheres.items():
            if (v.type == 'photosphere'):
                f = open('malla.grid', 'w')
                f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
                f.write("           b) The first six characters of the last line                \n")
                f.write("          in the header (if any) must contain the symbol ---       \n")
                f.write("\n")                                                                       
                f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
                f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
                f.write("-----------------------------------------------------------------------\n")

                ind_low = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[0])).argmin()
                ind_top = (np.abs(v.spectrum.wavelength_axis - v.wvl_range_lambda[1])).argmin()

                low = v.spectrum.wavelength_axis[ind_low]
                top = v.spectrum.wavelength_axis[ind_top] + 1e-3
                delta = (v.spectrum.wavelength_axis[1] - v.spectrum.wavelength_axis[0])

                ff = open('LINEAS', 'r')
                flines = ff.readlines()
                ff.close()

                for i in range(len(v.lines)):
                    for l in flines:
                        tmp = l.split()
                        index = int(tmp[0].split('=')[0])
                        if (index == v.lines[0]):
                            wvl = float(tmp[2])                    
                                            
                f.write("{0}            :  {1}, {2}, {3}\n".format(str(v.lines)[1:-1], 1e3*(low-wvl), 1e3*delta, 1e3*(top-wvl)))
                f.close()
                
                v.n_lambda = sir_code.init(v.index)

    def add_topology(self, atmosphere_order):
        """
        Add a new topology

        Parameters
        ----------
        topology : str
            Topology
        
        Returns
        -------
        None

        """

        # Transform the order to a list of lists        
        if (self.verbose):
            print('  - {0}'.format(atmosphere_order))

        vertical_order = atmosphere_order.split('->')        
        order = []
        for k in vertical_order:
            name = k.strip().replace('(','').replace(')','').split('+')
            name = [k.strip() for k in name]
            
            tmp = []
            for n in name:
                # tmp.append(self.atmospheres[n])
                tmp.append(n)
                self.atmospheres[n].active = True

            order.append(tmp)            

        self.order_atmospheres.append(order)
                
    def initialize_ff(self):
        """
        Normalize all filling factors so that they add to one to avoid later problems
        """

        for atmospheres in self.order_atmospheres:
            for order in atmospheres:
                total_ff = 0.0
                for atm in order:
                    total_ff += self.atmospheres[atm].ff

                for atm in order:                    
                    self.atmospheres[atm].ff /= total_ff                


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

    def synthesize_spectral_region(self, spectral_region, perturbation=False):
        """
        Synthesize all atmospheres

        Parameters
        ----------
        spectral_region : str
            Spectral region to synthesize
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """
        stokes = None
        stokes_out = None

        for atmospheres in self.order_atmospheres:
            
            for n, order in enumerate(atmospheres):
                                                                
                for k, atm in enumerate(order):
                                        
                    if (self.atmospheres[atm].spectrum.name == spectral_region):
                                                                                
                        if (n > 0):
                            ind_low, ind_top = self.atmospheres[atm].wvl_range
                            if (perturbation):
                                stokes_out = self.atmospheres[atm].spectrum.stokes_perturbed[:, ind_low:ind_top]
                            else:
                                stokes_out = self.atmospheres[atm].spectrum.stokes[:, ind_low:ind_top]
                                                                                
                        if (k == 0):
                            stokes = self.atmospheres[atm].ff * self.atmospheres[atm].synthesize(stokes_out)
                        else:                        
                            stokes += self.atmospheres[atm].ff * self.atmospheres[atm].synthesize(stokes_out)
                                    
                        ind_low, ind_top = self.atmospheres[atm].wvl_range
                        
                        if (perturbation):
                            self.atmospheres[atm].spectrum.stokes_perturbed[:, ind_low:ind_top+1] = stokes
                        else:
                            self.atmospheres[atm].spectrum.stokes[:, ind_low:ind_top+1] = stokes

    def synthesize(self, perturbation=False):
        """
        Synthesize all atmospheres

        Parameters
        ----------
        perturbation : bool
            Set to True if you are synthesizing with a perturbation. In this case, the synthesis
            is saved in spectrum.stokes_perturbed instead of spectrum.stokes
        
        Returns
        -------
        None

        """

        for k, v in self.spectrum.items():
            self.synthesize_spectral_region(k, perturbation=perturbation)


    def find_active_parameters(self, cycle):
        """
        Find all active parameters in all active atmospheres in the current cycle

        Parameters
        ----------
        cycle : int
            Cycle to consider
        
        Returns
        -------
        None

        """
        pars = []
        self.nodes = []
        left = 0
        right = 0
        for atmospheres in self.order_atmospheres:
            for n, order in enumerate(atmospheres):
                for k, atm in enumerate(order):                    
                    for l, par in self.atmospheres[atm].cycles.items():
                        #print(l, par)
                        if (par is not None):
                            if (hazel.util.isint(par[cycle])):
                                if (par[cycle] > 0):
                                    # [Atmosphere name, n_nodes, nodes, value, range]
                                    self.atmospheres[atm].nodes[l] = np.zeros(par[cycle])
                                    self.atmospheres[atm].n_nodes[l] = par[cycle]

                                    right += par[cycle]
                                    
                                    n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                    tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 
                                        'ranges': self.atmospheres[atm].ranges[l], 'delta': self.atmospheres[atm].epsilon[l],
                                        'left': left, 'right': right}

                                    self.nodes.append(self.atmospheres[atm].nodes[l])
                                        
                                    left = copy.copy(right)
                                                                                                            
                                    pars.append(tmp)
                            else:
                                self.atmospheres[atm].nodes[l] = np.zeros(par[cycle])
                                self.atmospheres[atm].n_nodes[l] = par[cycle]

                                right += par[cycle]

                                n_lambda = len(self.atmospheres[atm].spectrum.wavelength_axis)
                                tmp = {'atm': atm, 'n_nodes': par[cycle], 'parameter': l, 
                                    'ranges': self.atmospheres[atm].ranges[l], 'delta': self.atmospheres[atm].epsilon[l], 
                                    'left': left, 'right': right}

                                self.nodes.append(self.atmospheres[atm].nodes[l])

                                left = copy.copy(right)

                                pars.append(tmp)

        self.active_meta = pars        
        self.nodes = np.concatenate(self.nodes).ravel()
        
    def synthesize_and_compute_rf(self, compute_rf=False):
        """
        Compute response functions for all free parameters according to all active_parameters

        Parameters
        ----------
        cycle : int
            Cycle to consider
        
        Returns
        -------
        None

        """

        self.synthesize()

        if (not compute_rf):
            return

        n_active_pars = len(self.active_meta)

        loop = 0

        for par in self.active_meta:
            nodes = self.nodes[par['left']:par['right']]

            if (self.verbose):
                print(" * RF to {0} - {1} - nodes={2}".format(par['atm'], par['parameter'], par['n_nodes']))
                        
            for i in range(par['n_nodes']):
                perturbation = np.zeros(par['n_nodes'])
                if (nodes[i] == 0):
                    perturbation[i] = self.epsilon * par['delta']
                else:
                    perturbation[i] = self.epsilon * nodes[i]

                # Perturb this parameter
                self.atmospheres[par['atm']].nodes[par['parameter']] = nodes + perturbation
                
                # Synthesize
                self.synthesize(perturbation=True)

                # And come back to the original value of the nodes
                self.atmospheres[par['atm']].nodes[par['parameter']] = nodes

                rf = np.expand_dims((self.spectrum['spec1'].stokes - self.spectrum['spec1'].stokes_perturbed) / perturbation[i], 0)
                
                if (loop == 0):
                    self.response = rf
                else:
                    self.response = np.vstack([self.response, rf])

                loop += 1

    def flatten_parameters_to_reference(self):
        """
        Flatten all current parameters to the reference atmosphere

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """
        for k, v in self.atmospheres.items():
            v.set_reference_to_current_parameters()

    def set_new_model(self, nodes):
        """
        Set the nodes of the current model to the values passed on the arguments

        Parameters
        ----------
        nodes : float
            Array with the new set of nodes
        
        Returns
        -------
        None

        """

        n_active_pars = len(self.active_meta)
                
        for par in self.active_meta:
            left = par['left']
            right = par['right']

            self.atmospheres[par['atm']].nodes[par['parameter']] = nodes[left:right]

    def modified_svd_inverse(self, H, tol=1e-6):
        """
        Compute the inverse of the Hessian matrix using a modified SVD, by thresholding each subpsace separately

        Parameters
        ----------
        H : float
            Hessian matrix

        tol : float
            Tolerance for the singular value of each subspace
        
        Returns
        -------
        None

        """

        U, w, VT = np.linalg.svd(H, full_matrices=False)

        w_new = np.zeros_like(w)
        
        for par in self.active_meta:
            left = par['left']
            right = par['right']

            Ui = np.zeros_like(U)
            Ui[:,left:right] = U[:,left:right]

            Gamma_i = np.diagonal(np.diag(w) @ Ui.T @ U).copy()
            
            wmax = np.max(np.abs(Gamma_i))            
            Gamma_i[np.abs(Gamma_i) < tol*wmax] = 0.0

            w_new += Gamma_i

        w_new_inv = np.zeros_like(w)
        ind = np.where(w_new != 0)[0]
        w_new_inv[ind] = 1.0 / w_new[ind]
        
        return U, w_new_inv, VT
        

    def compute_chi2(self, obs, noise, only_chi2=False):
        """
        Compute chi2 for all spectral regions

        Parameters
        ----------
        obs : float
            Vector of observations
        only_chi2 : bool
            Control whether the gradient and Hessian is returned
        
        Returns
        -------
        None

        """
        for k, v in self.spectrum.items():
            residual = (v.stokes - obs)
            chi2 = np.sum(residual**2 * self.factor_chi2)            
            
            if (not only_chi2):
                dchi2 = -2.0 * np.sum(self.response * residual[None,:,:] * self.factor_chi2[None,:,:], axis=(1,2))            
                ddchi2 = 2.0 * np.sum(self.response[None,:,:,:] * self.response[:,None,:,:] * self.factor_chi2[None,None,:,:], axis=(2,3))
                return chi2, dchi2, ddchi2
            else:
                return chi2

    def invert(self):
        """
        Invert all atmospheres

        Parameters
        ----------
        None
        
        Returns
        -------
        None

        """

        self.spectrum['spec1'].observed_handle.open()
        obs, noise = self.spectrum['spec1'].observed_handle.read()

        dof = np.prod(obs.shape)

        self.factor_chi2 = 1.0 / (dof * noise**2)

        self.synthesize_and_compute_rf()
    
        f, ax = pl.subplots(nrows=2, ncols=2)
        ax = ax.flatten()
        for i in range(4):
            ax[i].plot(obs[i,:], '--')
            ax[i].plot(self.spectrum['spec1'].stokes[i,:], '.')

        lambdaLM = 10.0
        lambda_opt = 10.0
    
        for cycle in range(self.n_cycles):
            self.find_active_parameters(cycle)

            for iteration in range(5):
                
                self.synthesize_and_compute_rf(compute_rf=True)
                chi2, dchi2, ddchi2 = self.compute_chi2(obs, noise)

                if (self.verbose):
                    print('It: {0} - chi2: {1} - lambda: {2}'.format(iteration, chi2, lambda_opt))
                
                chi2_arr = []
                lambdas = []
                sols = []
                keepon = True
                loop = 0

                # Backtracking
                while keepon:                    

                    H = 0.5 * ddchi2
                    H += np.diag(lambdaLM * np.diag(H))
                    gradF = 0.5 * dchi2

                    U, w_inv, VT = self.modified_svd_inverse(H, tol=1e-8)

                    # xnew = xold - H^-1 * grad F
                    delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
                    
                    new_solution = self.nodes + delta
                    sols.append(new_solution)

                    self.set_new_model(new_solution)

                    self.synthesize_and_compute_rf()

                    chi2_arr.append(self.compute_chi2(obs, noise, only_chi2=True))
                    lambdas.append(lambdaLM)

                    if (self.verbose):
                        print('  - Backtracking: {0:2d} - lambda: {1:7.5f} - chi2: {2:7.5f}'.format(loop, lambdaLM, chi2_arr[-1]))
                    
                    if (loop > 0):
                        if (chi2_arr[-1] > chi2_arr[-2]):
                            keepon = False

                    if (lambdaLM < 1e-3):
                        keepon = False

                    lambdaLM /= np.sqrt(10.0)
                    loop += 1

                # Parabolic interpolation of the optimal value of lambda
                coeff = np.polyfit(np.log(lambdas[-3:]), chi2_arr[-3:], 2)
                lambda_opt = np.exp(-coeff[1] / (2.0*coeff[0]))

                if (self.verbose):
                    print('  * Optimal lambda: {0}'.format(lambda_opt))

                # Give the final step
                H = 0.5 * ddchi2
                H += np.diag(lambda_opt * np.diag(H))
                gradF = 0.5 * dchi2

                U, w_inv, VT = self.modified_svd_inverse(H, tol=1e-5)

                # xnew = xold - H^-1 * grad F
                delta = -VT.T.dot(np.diag(w_inv)).dot(U.T).dot(gradF)
                
                self.nodes += delta

                for i in range(4):
                    ax[i].plot(self.spectrum['spec1'].stokes[i,:])

                # Increase the optimal by 100 to find again the optimal value
                lambdaLM *= 100.0

            self.flatten_parameters_to_reference()
            stop()

    def plot_stokes(self):        
        for atmospheres in self.order_atmospheres:
            f, ax = pl.subplots(nrows=2, ncols=2)
            ax = ax.flatten()
            for i in range(4):
                ax[i].plot(atmospheres[0][-1].spectrum.stokes[i,:])
