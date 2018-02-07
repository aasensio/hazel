from collections import OrderedDict
import numpy as np
import os
from .util import i0_allen
from . import pyhazel
from . import pysir
from .hsra import hsra_continuum
from ipdb import set_trace as stop

__all__ = ['Hazel_atmosphere', 'SIR_atmosphere']
    
sir_parameters = OrderedDict.fromkeys('T B thetaB phiB v')

class General_atmosphere(object):
    def __init__(self, atm_type):
        self.ff = 1.0

        self.active_lines = []
        self.wavelength = dict()
        self.wavelength_range = dict()
        self.wvl_axis = dict()
        self.wvl_range = dict()
        self.type = atm_type
        self.spectrum = dict()

        self.multiplets = {'10830': 10829.0911, '3888': 3888.6046, '7065': 7065.7085, '5876': 5875.9663}

    # def set_

class Parametric_atmosphere(General_atmosphere):
    def __init__(self):
    
        super().__init__('parametric')

        pass


class Hazel_atmosphere(General_atmosphere):
    def __init__(self):
    
        super().__init__('chromosphere')

        self.parameters = dict()
        

    # def add_active_line(self, line, spectrum, wvl_range):
    #     """
    #     Add an active lines in this atmosphere
        
    #     Parameters
    #     ----------
    #     lines : str
    #         Line to activate: ['10830','5876']
    #     spectrum : Spectrum
    #         Spectrum object
    #     wvl_range : float
    #         Vector containing wavelength range over which to synthesize this line
        
    #     Returns
    #     -------
    #     None
    
    #     """        
    #     self.active_lines.append(line)
    #     self.wavelength_range[line] = wvl_range
    #     ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
    #     ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

    #     self.spectrum[line] = spectrum
    #     self.wvl_axis[line] = spectrum.wavelength_axis[ind_low:ind_top]
    #     self.wvl_range[line] = [ind_low, ind_top]

    #     self.parameters[line] = OrderedDict()
    #     self.parameters[line]['B'] = 0.0
    #     self.parameters[line]['thetaB'] = 0.0
    #     self.parameters[line]['phiB'] = 0.0
    #     self.parameters[line]['h'] = 3.0
    #     self.parameters[line]['tau'] = 1.0
    #     self.parameters[line]['v'] = 0.0
    #     self.parameters[line]['deltav'] = 8.0
    #     self.parameters[line]['beta'] = 1.0
    #     self.parameters[line]['a'] = 0.0

    def add_active_line(self, line, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------
        lines : str
            Line to activate: ['10830','5876']
        spectrum : Spectrum
            Spectrum object
        wvl_range : float
            Vector containing wavelength range over which to synthesize this line
        
        Returns
        -------
        None
    
        """        
        self.active_line = line
        self.wavelength_range = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top]

        self.parameters = OrderedDict()
        self.parameters['B'] = 0.0
        self.parameters['thetaB'] = 0.0
        self.parameters['phiB'] = 0.0
        self.parameters['h'] = 3.0
        self.parameters['tau'] = 1.0
        self.parameters['v'] = 0.0
        self.parameters['deltav'] = 8.0
        self.parameters['beta'] = 1.0
        self.parameters['a'] = 0.0

    def synthesize(self, stokes=None):
        
        synModeInput = 5
        nSlabsInput = 1
        B1Input = np.asarray([self.parameters['B'], self.parameters['thetaB'], self.parameters['phiB']])
        B2Input = np.asarray([0.0,0.0,0.0])
        hInput = self.parameters['h']
        tau1Input = self.parameters['tau']
        tau2Input = 0.e0        
        transInput = 1
        atomicPolInput = 1
        magoptInput = 1
        anglesInput = np.asarray([0.0,0.0,90.0])
        lambdaAxisInput = self.wvl_axis - self.multiplets[self.active_line]
        nLambdaInput = len(lambdaAxisInput)
        
        if (stokes is None):
            boundaryInput  = np.zeros((nLambdaInput,4))
            boundaryInput[:,0] = hsra_continuum(self.multiplets[self.active_line]) #i0_allen(self.multiplets[self.active_line],1.0)            
        else:            
            boundaryInput = stokes.T * hsra_continuum(self.multiplets[self.active_line]) #i0_allen(self.multiplets[self.active_line],1.0)
            
        dopplerWidthInput = self.parameters['deltav']
        dopplerWidth2Input = 0.e0
        dampingInput = self.parameters['a']
        dopplerVelocityInput = self.parameters['v']
        dopplerVelocity2Input = 0.e0
        ffInput = 0.e0
        betaInput = 1.0
        beta2Input = 1.0
        nbarInput = np.asarray([0.0,0.0,0.0,0.0])
        omegaInput = np.asarray([0.0,0.0,0.0,0.0])
        normalization = 0
        
        args = (synModeInput, nSlabsInput, 
            B1Input, B2Input, hInput, tau1Input, tau2Input, boundaryInput, transInput, 
            atomicPolInput, magoptInput, anglesInput, nLambdaInput, lambdaAxisInput, dopplerWidthInput, 
            dopplerWidth2Input, dampingInput, dopplerVelocityInput, dopplerVelocity2Input, 
            ffInput, betaInput, beta2Input, nbarInput, omegaInput, normalization)

        l, stokes, eta, eps = pyhazel._synth(*args)
        
        return stokes / hsra_continuum(self.multiplets[self.active_line])

class SIR_atmosphere(General_atmosphere):
    def __init__(self):
        
        super().__init__('photosphere')
        self.ff = 1.0
        self.initialization = True
        self.parameters = dict()
        
    def list_lines(self):
        """
        List the lines available in SIR for synthesis
            
        """
        f = open('LINEAS', 'r')
        lines = f.readlines()
        f.close()

        print("Available lines:")
        for l in lines[:-1]:
            print(l[:-1])

    def add_active_line(self, lines, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------
        lines : str
            Line to activate: ['10830','5876']
        spectrum : Spectrum
            Spectrum object
        wvl_range : float
            Vector containing wavelength range over which to synthesize this line
        
        Returns
        -------
        None
    
        """

        if (self.initialization):

            f = open('malla.grid', 'w')
            f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
            f.write("           b) The first six characters of the last line                \n")
            f.write("          in the header (if any) must contain the symbol ---       \n")
            f.write("\n")                                                                       
            f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
            f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
            f.write("-----------------------------------------------------------------------\n")

            self.initialization = False
        else:
            f = open('malla.grid', 'a')


        # self.wavelength_range[line] = wvl_range
        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = [ind_low, ind_top]

        low = spectrum.wavelength_axis[ind_low]
        top = spectrum.wavelength_axis[ind_top] + 1e-3
        delta = (spectrum.wavelength_axis[1] - spectrum.wavelength_axis[0])

        ff = open('LINEAS', 'r')
        flines = ff.readlines()
        ff.close()

        for i in range(len(lines)):
            for l in flines:
                tmp = l.split()
                index = int(tmp[0].split('=')[0])
                if (index == lines[0]):
                    wvl = float(tmp[2])                    
                                    
        f.write("{0}            :  {1}, {2}, {3}\n".format(str(lines)[1:-1], 1e3*(low-wvl), 1e3*delta, 1e3*(top-wvl)))
        f.close()


        # if (not os.path.exists('LINEAS')):
        #     local = str(__file__).split('/')
        #     sdir = '/'.join(local[0:-2])+'/data'
        #     shutil.copy(sdir+'/LINEAS', os.getcwd())

        # if (not os.path.exists('THEVENIN')):
        #     local = str(__file__).split('/')
        #     sdir = '/'.join(local[0:-2])+'/data'
        #     shutil.copy(sdir+'/THEVENIN', os.getcwd())
            
        pysir.init()
            
    def _interpolate_nodes(logTau, nodes, nodes_logtau=None, variable=None):
        """
        Generate a model atmosphere by interpolating the defined nodes. The interpolation
        order depends on the number of nodes.
        
        Parameters
        ----------
        logTau : float
            Vector of log optical depth at 500 nm
        nodes : float
            List with the position of the nodes

        Returns
        -------
        real
            Vector with the interpolated atmosphere
    
        """
        n = logTau.shape[0]

        if (variable is None):
            variable = np.zeros_like(logTau)

        out = np.zeros_like(logTau)

        if (nodes_logtau is not None):
            ind = np.argsort(nodes_logtau)
            spl = PchipInterpolator(nodes_logtau[ind], nodes[ind], extrapolate=True)
            out = spl(logTau)
            return out

        if (len(nodes) == 1):
            out = variable + nodes
        
        if (len(nodes) >= 2):
            pos = np.linspace(0, n-1, len(nodes), dtype=int)
            coeff = np.polyfit(logTau[pos], nodes, len(nodes)-1)
            out = variable + np.polyval(coeff, logTau)
        return out

    def load_model(self, model_file):
        out = np.loadtxt('model.mod', dtype=np.float32)
        out = np.delete(out, 2, axis=1)    
        self.set_model(out)

    def set_model(self, model, macroturbulence=0.0, fillingFactor=1.0, stray=0.0):
        if (model.shape[1] == 7):
            model = np.insert(model, 2, -np.ones(model.shape[0]), axis=1)
    # Boundary condition for Pe
            model[-1,2] = 1.11634e-01

        self.model = model
        self.macroturbulence = macroturbulence
        self.filling_factor = fillingFactor
        self.stray = stray

    def synthesize(self, returnRF=False):
        """Carry out the synthesis and returns the Stokes parameters and the response 
        functions to all physical variables at all depths
        
        Args:
            model (float array): an array of size [nDepth x 7] or [nDepth x 8], where the columns contain the depth stratification of:
                - log tau
                - Temperature [K]
                - Electron pressure [dyn cm^-2]  (optional)
                - Microturbulent velocity [km/s]
                - Magnetic field strength [G]
                - Line-of-sight velocity [km/s]
                - Magnetic field inclination [deg]
                - Magnetic field azimuth [deg]
            macroturbulence (float, optional): macroturbulence velocity [km/s]. Default: 0
            fillingFactor (float, optional): filling factor. Default: 1
            stray (float, optional): stray light in %. Default: 0
            returnRF (bool, optional): return response functions
        
        Returns:
            stokes: (float array) Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (5,nLambda)
            rf: (float array) Response functions to T, Pe, vmic, B, v, theta, phi, all of size (4,nLambda,nDepth), plus the RF to macroturbulence of size (4,nLambda)
                            It is not returned if returnRF=False
        """
        
        if (returnRF):
            stokes, rf = pysir.synthRF(self.model, self.macroturbulence, self.filling_factor, self.stray)
            return stokes[1:,:], rf
        else:
            stokes = pysir.synth(self.model, self.macroturbulence, self.filling_factor, self.stray)
            return stokes[1:,:]

        return stokes[1:,:], rf