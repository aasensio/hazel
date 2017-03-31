Configuration
=============

The code is controlled by a single human-readable configuration file. 
In order to use this option, you need to have the ``configparser``
package installed in your system. It can be downloaded from
``http://www.voidspace.org.uk/python/configobj.html`` or installed as a standard
Python package using ``pip install configparser''. The serial version of the code
is run using:

::

    ./run.py conf.ini

and the parallel code is run with

::

    ./run.py conf.ini nProcessors

An example of the file, that is self-explanatory, is:

::

    # Hazel configuration File

    #####################
    # General information
    #####################

    [Files]
    Input model file = 'ATOMS/helium.mod'
    File with observations = 'OBSERVATION/test_2comp.prof'
    File with inverted profiles = 'test.inversion'
    File with inverted parameters = 'test.parameters'

    [Working mode]
    Action = 'inversion'                        # 'synthesis' or 'inversion'
    Verbose = no
    Linear system solver = 'LU'                 # 'LU' or 'CG'
    Stopping volume for DIRECT = 0.001

    [General parameters]
    Synthesis mode = 'exact'                    # 'thin' or 'exact'
    Include stimulated emission = yes
    Include magnetic field  = yes
    Include Paschen-Back effect = yes
    Include atomic level polarization = yes
    Include magneto-optical effects in the RT = yes
    Include stimulated emission in the RT = yes
    Multiplet = 10830                           # 10830, 5876, 7065, 3889 A
    Line-of-sight angles = 0.0, 0.0, 90.0       # theta, chi, gamma deg
    Wavelength axis = -3.0, 2.5, 200            # Minimum, maximum and number of grid points

    #####################
    # Synthesis parameters
    #####################
    [Synthesis]
    Number of slabs = '1'                       # '1' -> single slab, '1+1' -> two slabs with same field, '1+1B' -> 2 slabs with different field, '2' -> two slabs added with a filling factor
    Boundary condition = 4.098e-5, 0.0, 0.0, 0.0      # I0, Q0, U0, V0
    a = 0.0
    height = 3.0                                # Real height if positive, apparent height if negative arcsec
    ff = 0.0
        [[Slab 1]]
        B =         0.0         # G
        thetaB =    0.0         # deg
        chiB =      0.0         # deg
        vdopp =     8.0         # km/s
        tau =       1.0
        vmac =      0.0         # Positive is redshift km/s
        beta =      1.0
        [[Slab 2]]
        B =         0.0         # G
        thetaB =    0.0         # deg
        chiB =      0.0         # deg
        vdopp =     0.0         # km/s
        tau =       0.0
        vmac =      0.0         # Positive is redshift km/s
        beta =      1.0

    #####################
    # Ranges for the DIRECT method [min, max]
    #####################
    [Ranges]
    a =             0,0.5
    ff =            0.0,1.0
        [[Slab 1]]
        B =         800,1100
        thetaB =    0,180
        chiB =      0,180
        vdopp =     2,12
        tau =       0.1,2
        vmac =      -5,5
        beta =      0.5,2
        [[Slab 2]]
        B =         800,1100
        thetaB =    0,180
        chiB =      0,180
        vdopp =     2,12
        tau =       0.1,2
        vmac =      -5,5
        beta =      0.5,2
        
    #####################
    # Parameters to invert
    #####################
    [Inversion]
    Iterations in LM = 20
    Number of cycles = 4
    Inversion modes = 'DIRECT', 'LM', 'DIRECT', 'LM'        # 'DIRECT' for DIRECT algorithm and 'LM' for Levenberg-Marquardt
        [[Cycles]]
        a =             1, 1, 0, 0
        ff =            0, 0, 0, 0
            [[[Slab 1]]]
            B =         0, 0, 1, 1
            thetaB =    0, 0, 1, 1
            chiB =      0, 0, 1, 1
            vdopp =     1, 1, 0, 0
            tau =       1, 1, 0, 0
            vmac =      1, 1, 0, 0
            beta =      0, 0, 0, 0
            [[[Slab 2]]]
            B =         0, 0, 0, 0
            thetaB =    0, 0, 0, 0
            chiB =      0, 0, 0, 0
            vdopp =     0, 0, 0, 0
            tau =       0, 0, 0, 0
            vmac =      0, 0, 0, 0
            beta =      0, 0, 0, 0
        [[Weights]]
            Stokes I =  1.0, 1.0, 1.0, 1.0
            Stokes Q =  0.0, 0.0, 1.0, 1.0
            Stokes U =  0.0, 0.0, 1.0, 1.0
            Stokes V =  0.0, 0.0, 1.0, 1.0


The code will convolve with a spectral point spread function if a file called ``psf.txt`` is placed on the
directory where you run the code. This file contains a first line with the number of points of the PSF.
Then a single column giving the transmission in a wavelength axis with a step equal to that of the observations.
If this file is not present, no convolution will be done.