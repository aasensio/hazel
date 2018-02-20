from collections import OrderedDict
import numpy as np
import os
from hazel.util import i0_allen
from hazel.codes import hazel_code, sir_code
from hazel.hsra import hsra_continuum
from ipdb import set_trace as stop
from hazel.sir import Sir


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
        self.active = False

        self.multiplets = {'10830': 10829.0911, '3888': 3888.6046, '7065': 7065.7085, '5876': 5875.9663}


class Parametric_atmosphere(General_atmosphere):
    def __init__(self):
    
        super().__init__('parametric')

        pass


class Hazel_atmosphere(General_atmosphere):
    def __init__(self):
    
        super().__init__('chromosphere')

        self.parameters = dict()
        
    def add_active_line(self, index, line, spectrum, wvl_range):
        """
        Add an active lines in this atmosphere
        
        Parameters
        ----------
        index : int
            Index of this chromosphere
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
        self.wvl_range = [ind_low, ind_top+1]

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

        self.cycles = OrderedDict()
        self.cycles['B'] = None
        self.cycles['thetaB'] = None
        self.cycles['phiB'] = None
        self.cycles['h'] = None
        self.cycles['tau'] = None
        self.cycles['v'] = None
        self.cycles['deltav'] = None
        self.cycles['beta'] = None
        self.cycles['a'] = None

        self.index = index + 1

    def set_parameters(self, pars):
        self.parameters['B'] = pars[0]
        self.parameters['thetaB'] = pars[1]
        self.parameters['phiB'] = pars[2]        
        self.parameters['tau'] = pars[3]
        self.parameters['v'] = pars[4]
        self.parameters['deltav'] = pars[5]
        self.parameters['beta'] = pars[6]
        self.parameters['a'] = pars[7]

    def load_reference_model(self, model_file, verbose):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            verbosity flag

        Returns
        -------
            None
        """
        extension = os.path.splitext(model_file)[1][1:]
        if (extension == '1d'):
            if (verbose):
                print('  - Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            out = np.loadtxt(model_file, skiprows=1)            
            self.set_parameters(out)
        
        if (extension == 'h5'):
            if (verbose):
                print('  - Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            self.model_file = model_file

    def synthesize(self, stokes=None):
        
        B1Input = np.asarray([self.parameters['B'], self.parameters['thetaB'], self.parameters['phiB']])
        hInput = self.parameters['h']
        tau1Input = self.parameters['tau']
        transInput = 1
        anglesInput = np.asarray([0.0,0.0,90.0])
        lambdaAxisInput = self.wvl_axis - self.multiplets[self.active_line]
        nLambdaInput = len(lambdaAxisInput)
        
        if (stokes is None):
            boundaryInput  = np.asfortranarray(np.zeros((4,nLambdaInput)))
            boundaryInput[0,:] = hsra_continuum(self.multiplets[self.active_line]) #i0_allen(self.multiplets[self.active_line],1.0)            
        else:            
            boundaryInput = np.asfortranarray(stokes * hsra_continuum(self.multiplets[self.active_line])) #i0_allen(self.multiplets[self.active_line],1.0)
                    
        dopplerWidthInput = self.parameters['deltav']
        dampingInput = self.parameters['a']
        dopplerVelocityInput = self.parameters['v']
        betaInput = self.parameters['beta']
        nbarInput = np.asarray([0.0,0.0,0.0,0.0])
        omegaInput = np.asarray([0.0,0.0,0.0,0.0])
        
        args = (self.index, B1Input, hInput, tau1Input, boundaryInput, transInput, 
            anglesInput, nLambdaInput, lambdaAxisInput, dopplerWidthInput, 
            dampingInput, dopplerVelocityInput, 
            betaInput, nbarInput, omegaInput)
        
        l, stokes = hazel_code._synth(*args)
        
        return stokes / hsra_continuum(self.multiplets[self.active_line])

class SIR_atmosphere(General_atmosphere):
    def __init__(self):
        
        super().__init__('photosphere')
        self.ff = 1.0        
        self.parameters = dict()

        self.cycles = OrderedDict()
        self.cycles['T'] = None
        self.cycles['vmic'] = None
        self.cycles['v'] = None
        self.cycles['B'] = None
        self.cycles['thetaB'] = None
        self.cycles['phiB'] = None
        
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

    def add_active_line(self, index, lines, spectrum, wvl_range, append=False):
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
        
        # if (not append):

        f = open('malla.grid', 'w')
        f.write("IMPORTANT: a) All items must be separated by commas.                 \n")
        f.write("           b) The first six characters of the last line                \n")
        f.write("          in the header (if any) must contain the symbol ---       \n")
        f.write("\n")                                                                       
        f.write("Line and blends indices   :   Initial lambda     Step     Final lambda \n")
        f.write("(in this order)                    (mA)          (mA)         (mA)     \n")
        f.write("-----------------------------------------------------------------------\n")
            
        # else:
            # f = open('malla.grid', 'a')      

        
        self.index = index + 1

        ind_low = (np.abs(spectrum.wavelength_axis - wvl_range[0])).argmin()
        ind_top = (np.abs(spectrum.wavelength_axis - wvl_range[1])).argmin()

        self.spectrum = spectrum
        self.wvl_axis = spectrum.wavelength_axis[ind_low:ind_top+1]
        self.wvl_range = np.array([ind_low, ind_top+1])
        
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

        self.n_lambda = sir_code.init(self.index)
            
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

    def load_reference_model(self, model_file, verbose):
        """
        Load a reference model or a model for every pixel for synthesis/inversion

        Parameters
        ----------
        model_file : str
            String with the name of the file. Extensions can currently be "1d" or "h5"
        verbose : bool
            Verbosity

        Returns
        -------
        None
        """
        extension = os.path.splitext(model_file)[1][1:]
        if (extension == '1d'):
            if (verbose):
                print('  - Reading 1D model {0} as reference'.format(model_file))
            self.model_type = '1d'
            out = np.loadtxt(model_file, dtype=np.float32, skiprows=1)
            out = np.delete(out, 2, axis=1)
            self.set_model(out)
        
        if (extension == 'h5'):
            if (verbose):
                print('  - Reading 3D model {0} as reference'.format(model_file))
            self.model_type = '3d'
            self.model_file = model_file
            
    def set_parameters(self, model, macroturbulence=0.0, fillingFactor=1.0, stray=0.0):
        if (model.shape[1] == 7):
            model = np.insert(model, 2, -np.ones(model.shape[0]), axis=1)
    # Boundary condition for Pe
            model[-1,2] = 1.11634e-01

        self.model = model
        self.macroturbulence = macroturbulence
        self.filling_factor = fillingFactor
        self.stray = stray

    def synthesize(self, stokes_in, returnRF=False):
        """Carry out the synthesis and returns the Stokes parameters and the response 
        functions to all physical variables at all depths
        
        Parameters
        ----------
        model : float
            An array of size [nDepth x 7] or [nDepth x 8], where the columns contain the depth stratification of:
                - log tau
                - Temperature [K]
                - Electron pressure [dyn cm^-2]  (optional)
                - Microturbulent velocity [km/s]
                - Magnetic field strength [G]
                - Line-of-sight velocity [km/s]
                - Magnetic field inclination [deg]
                - Magnetic field azimuth [deg]
        macroturbulence : float, optional
            Macroturbulence velocity [km/s]. Default: 0
        fillingFactor : float, optional
            Filling factor. Default: 1
        stray : float, optional
            Stray light in %. Default: 0
        returnRF : bool, optional
            Return response functions
        
        Returns
        -------
        
        stokes : float
            Stokes parameters, with the first index containing the wavelength displacement and the remaining
                                    containing I, Q, U and V. Size (5,nLambda)
        rf: float 
            Response functions to T, Pe, vmic, B, v, theta, phi, all of size (4,nLambda,nDepth), plus the RF to macroturbulence of size (4,nLambda)
                            It is not returned if returnRF=False
        """
        if (returnRF):
            # stokes, rf = self.sir.synthRF(self.model, self.macroturbulence, self.filling_factor, self.stray)
            stokes, rf = sir_code.synthRF(self.index, self.n_lambda, self.model, self.macroturbulence, self.filling_factor, self.stray)
            return stokes[1:,:], rf
        else:
            # stokes = self.sir.synth(self.model, self.macroturbulence, self.filling_factor, self.stray)
            stokes = sir_code.synth(self.index, self.n_lambda, self.model, self.macroturbulence, self.filling_factor, self.stray)
            return stokes[1:,:]

        return stokes[1:,:], rf