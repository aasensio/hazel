#!/usr/bin/env python
from __future__ import print_function
from configobj import ConfigObj
import sys
import os
from subprocess import call

def lower_to_sep(string, separator='='):
	line=string.partition(separator)
	string=str(line[0]).lower()+str(line[1])+str(line[2])
	return string


if (len(sys.argv) < 2):
	print("Example usage: ./run.py conf.ini n_procs")
	exit()



print("Using configuration file = "+sys.argv[1])

# Transform all keys to lowercase to avoid problems with
# upper/lower case
f = open(sys.argv[1],'r')
input_lines = f.readlines()
f.close()
input_lower = ['']
for l in input_lines:
	input_lower.append(lower_to_sep(l)) # Convert keys to lowercase
        
# Parse configuration file
config = ConfigObj(input_lower)

id = ''
if (len(sys.argv) == 3):
	id = sys.argv[2] 

#*********************************
# Save config_inversion.dat file
#*********************************
f = open('config_inversion.dat','w')
f.write('***********************************************\n')
f.write('* Configuration file for the multi-term program\n')
f.write('***********************************************\n')
f.write('\n')
f.write("# Input model file\n")
f.write("'"+config['files']['input model file']+"'\n")
f.write('\n')
f.write('# Initial parameters file\n')
f.write("'init_parameters.dat'\n")
f.write('\n')
f.write('# Range of parameters for the DIRECT method\n')
f.write("'direct_range.dat'\n")
f.write('\n')
f.write("# Output file for the upper level rho^K_Q(J,J') in the reference frame of the vertical\n")
f.write("'ATOMIC_POL/vertical_upper.rho'\n")
f.write('\n')
f.write("# Output file for the lower level rho^K_Q(J,J') in the reference frame of the vertical\n")
f.write("'ATOMIC_POL/vertical_lower.rho'\n")
f.write('\n')
f.write("# Output file for the upper level rho^K_Q(J,J') in the reference frame of the magnetic field\n")
f.write("'ATOMIC_POL/magnetic_upper.rho'\n")
f.write('\n')
f.write("# Output file for the lower level rho^K_Q(J,J') in the reference frame of the magnetic field\n")
f.write("'ATOMIC_POL/magnetic_lower.rho'\n")
f.write('\n')
f.write("# Output absorption/emission coefficients\n")
f.write("'INVERTED/rtcoef.emer'\n")
f.write('\n')
f.write("# Output absorption/emission coefficients neglecting atomic polarization\n")
f.write("'INVERTED/rtcoef_noatompol.emer'\n")
f.write('\n')
f.write("# File with the observed profiles\n")
f.write("'"+config['files']['file with observations']+"'\n")
f.write('\n')
f.write("# File with the inverted profiles\n")
f.write("'"+config['files']['file with inverted profiles']+"'\n")
f.write('\n')
f.write("# File with the parameters from the inversion\n")
f.write("'"+config['files']['file with inverted parameters']+"'\n")
f.write('\n')
f.write("# File with the final errors\n")
f.write("'"+config['files']['file with errors in inverted parameters']+"'\n")
f.write('\n')
f.write("# File that sets the parameters to invert\n")
f.write("'invert_parameters.dat'\n")
f.write('\n')
f.write("# Verbose mode (0-> no, 1-> yes)\n")
if (config['working mode']['verbose'] == 'yes'):
	f.write('1\n')
else:
	f.write('0\n')
f.write('\n')
f.write("# Linear system solver (0-> LU, 1-> CG)\n")
if (config['working mode']['linear system solver'] == 'LU'):
	f.write('0\n')
else:
	f.write('1\n')
f.write('\n')
f.write("# Optically thin (0), slab no-MO (1), M-E (2), slab DELOPAR (3), simplified slab (4), exact slab (5)\n")
if (config['general parameters']['synthesis mode'] == 'exact'):
	f.write('5\n')
elif (config['general parameters']['synthesis mode'] == 'thin'):
	f.write('0\n')
else:
	print("Synthesis mode not supported : {0}".format(config['general parameters']['synthesis mode']))
	sys.exit()
f.write('\n')
f.write("# Synthesis mode -> 0 , Inversion mode -> 1\n")
if (config['working mode']['action'] == 'synthesis'):
	f.write('0')
elif (config['working mode']['action'] == 'inversion'):
	f.write('1')
else:
	print("Action mode not supported : {0}".format(config['working mode']['action']))
	sys.exit()
f.close()


#*********************************
# Save direct_range.dat file
#*********************************
f = open('direct_range.dat','w')
f.write("***********************************************\n")
f.write("* Ranges for the DIRECT method\n")
f.write("***********************************************\n")
f.write("\n")
f.write("# Maximum number of function evaluations (<0 -> don't use this criteria)\n")
f.write("-1\n")
f.write("\n")
f.write("# Reduction in the volume (<0 -> don't use this criteria, typically 0.01)\n")
f.write(config['working mode']['stopping volume for direct']+"\n")
f.write("\n")
f.write("# Magnetic field (0-Bmax)\n")
f.write("  ".join(config['ranges']['slab 1']['b'])+"\n")
f.write("\n")
f.write("# thetab  (0 .. 180)\n")
f.write("  ".join(config['ranges']['slab 1']['thetab'])+"\n")
f.write("\n")
f.write("# chib (0 .. 180)\n")
f.write("  ".join(config['ranges']['slab 1']['chib'])+"\n")
f.write("\n")
f.write("# vdopp (0 .. 20)\n")
f.write("  ".join(config['ranges']['slab 1']['vdopp'])+"\n")
f.write("\n")
f.write("# dtau (0 .. 5)\n")
f.write("  ".join(config['ranges']['slab 1']['tau'])+"\n")
f.write("\n")
f.write("# delta_collision (0 .. 18)\n")
f.write("0.d0  18.d0\n")
f.write("\n")
f.write("# vmacro (-10 .. 10)\n")
f.write("  ".join(config['ranges']['slab 1']['vmac'])+"\n")
f.write("\n")
f.write("# damping (0 .. 4)\n")
f.write("  ".join(config['ranges']['a'])+"\n")
f.write("\n")
f.write("# beta (0 .. 10)\n")
f.write("  ".join(config['ranges']['slab 1']['beta'])+"\n")
f.write("\n")
f.write("# height (0 .. 100)\n")
f.write("0.d0  100.d0\n")
f.write("\n")
f.write("# dtau2 (0 .. 5)\n")
f.write("  ".join(config['ranges']['slab 2']['tau'])+"\n")
f.write("\n")
f.write("# vmacro2 (-10 .. 10)\n")
f.write("  ".join(config['ranges']['slab 2']['vmac'])+"\n")
f.write("\n")
f.write("# Magnetic field 2 (0-Bmax)\n")
f.write("  ".join(config['ranges']['slab 2']['b'])+"\n")
f.write("\n")
f.write("# thetab 2 (0 .. 180)\n")
f.write("  ".join(config['ranges']['slab 2']['thetab'])+"\n")
f.write("\n")
f.write("# chib 2 (0 .. 180)\n")
f.write("  ".join(config['ranges']['slab 2']['chib'])+"\n")
f.write("\n")
f.write("# vdopp 2 (0 .. 20)\n")
f.write("  ".join(config['ranges']['slab 2']['vdopp'])+"\n")
f.write("\n")
f.write("# ff\n")
f.write("  ".join(config['ranges']['ff'])+"\n")
f.write("\n")
f.write("# beta2 (0 .. 10)\n")
f.write("  ".join(config['ranges']['slab 2']['beta'])+"\n")
f.close()

#*********************************
# Save invert_parameters.dat file
#*********************************
method = {'DIRECT': '2', 'LM': '1'}
f = open('invert_parameters.dat','w')
f.write("***********************************************\n")
f.write("* File defining the parameters to invert\n")
f.write("***********************************************\n")
f.write("\n")
f.write("# Maximum number of iterations\n")
f.write(config['inversion']['iterations in lm']+"\n")
f.write("\n")
f.write("# Number of cycles\n")
f.write(config['inversion']['number of cycles']+"\n")
f.write("\n")
f.write("# Invert the magnetic field strength\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['b'])+"\n")
f.write("\n")
f.write("# Invert the magnetic field inclination\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['thetab'])+"\n")
f.write("\n")
f.write("# Invert the magnetic field azimuth\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['chib'])+"\n")
f.write("\n")
f.write("# Invert the Doppler width\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['vdopp'])+"\n")
f.write("\n")
f.write("# Invert the optical depth or strength of the line\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['tau'])+"\n")
f.write("\n")
f.write("# Invert the D^2 of the lower level\n")
f.write(int(config['inversion']['number of cycles'])*" 0" + "\n")
f.write("\n")
f.write("# Invert the macroscopic velocity\n")
f.write(" ".join(config['inversion']['cycles']['slab 1']['vmac'])+"\n")
f.write("\n")
f.write("# Invert the damping\n")
f.write(" ".join(config['inversion']['cycles']['a'])+"\n")
f.write("\n")
f.write("# Invert beta\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['beta'])+"\n")
f.write("\n")
f.write("# Invert the height of the He atoms\n")
f.write(int(config['inversion']['number of cycles'])*" 0" + "\n")
f.write("\n")
f.write("# Invert the optical depth or strength of the line of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['tau'])+"\n")
f.write("\n")
f.write("# Invert the macroscopic velocity of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['vmac'])+"\n")
f.write("\n")
f.write("# Invert the magnetic field strength of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['b'])+"\n")
f.write("\n")
f.write("# Invert the magnetic field inclination of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['thetab'])+"\n")
f.write("\n")
f.write("# Invert the magnetic field azimuth of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['chib'])+"\n")
f.write("\n")
f.write("# Invert the Doppler width of component 2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['vdopp'])+"\n")
f.write("\n")
f.write("# Invert filling factor\n")
f.write(" ".join(config['inversion']['cycles']['ff'])+"\n")
f.write("\n")
f.write("# Invert beta2\n")
f.write(" ".join(config['inversion']['cycles']['slab 2']['beta'])+"\n")
f.write("\n")
f.write("# Weights for Stokes I in each cycle\n")
f.write(" ".join(config['inversion']['weights']['stokes i'])+"\n")
f.write("\n")
f.write("# Weights for Stokes Q in each cycle\n")
f.write(" ".join(config['inversion']['weights']['stokes q'])+"\n")
f.write("\n")
f.write("# Weights for Stokes U in each cycle\n")
f.write(" ".join(config['inversion']['weights']['stokes u'])+"\n")
f.write("\n")
f.write("# Weights for Stokes V in each cycle\n")
f.write(" ".join(config['inversion']['weights']['stokes v'])+"\n")
f.write("\n")
f.write("# Inversion modes (1-> LM, 2-> DIRECT, 3-> PIKAIA)\n")
f.write(" ".join([method[i] for i in config['inversion']['inversion modes']]))
f.close()

#*********************************
# Save init_parameters.dat file
#*********************************

# Define some dictionaries to facilitate writing the file
yesno = {'yes': '1', 'no': '0'}
slabs = {'1': '1', '1+1': '2', '1+1B': '3', '2': '-2'}
multiplets = {'10830': '1', '3889': '2', '7065': '3', '5876': '4'}
lambda0 = {'10830': '10829.0911', '3889': '3888.6046', '7065': '7065.7085', '5876': '5876.9663'}
if (config['synthesis']['number of slabs'] == '1'):
	nSlabs = 1
else:
	nSlabs = 2
includeff = False
if (config['synthesis']['number of slabs'] == '2'):
	includeff = True
	
f = open('init_parameters.dat','w')
f.write("***********************************************\n")
f.write("* File defining the specific experiment to solve\n")
f.write("***********************************************\n")
f.write("\n")
f.write("# Include stimulated emission (0-> no, 1-> yes)\n")
f.write(yesno[config['general parameters']['include stimulated emission']]+"\n")
f.write("\n")
f.write("# Include magnetic field (0-> no, 1-> yes)\n")
f.write(yesno[config['general parameters']['include magnetic field']]+"\n")
f.write("\n")
f.write("# Include depolarization rates (0-> no, 1-> yes)\n")
f.write("0\n")
f.write("\n")
f.write("# Value of delta if depolarization rates are included (not used if the previous value is 0)\n")
f.write("0.0\n")
f.write("\n")
f.write("# Include Paschen-Back effect (0-> no, 1-> yes)\n")
f.write(yesno[config['general parameters']['include paschen-back effect']]+"\n")
f.write("\n")
f.write("# Number of slabs (1-> 1 slab, 2-> 2 slabs with same B, 3-> 2 slabs with different B (add field below))\n")
f.write(slabs[config['synthesis']['number of slabs']]+"\n")
f.write("\n")
f.write("# Magnetic field strength [G], thetaB [degrees], chiB [degrees]\n")
if (nSlabs == 1):
	f.write(config['synthesis']['slab 1']['b']+"  "+config['synthesis']['slab 1']['thetab']+"  "+config['synthesis']['slab 1']['chib']+"\n")
else:
	f.write(config['synthesis']['slab 1']['b']+"  "+config['synthesis']['slab 1']['thetab']+"  "+config['synthesis']['slab 1']['chib']+"   "+\
		config['synthesis']['slab 2']['b']+"  "+config['synthesis']['slab 2']['thetab']+"  "+config['synthesis']['slab 2']['chib']+"\n")
f.write("\n")
f.write("# Apparent height (if negative) or real height (if positive) of the atoms in arcsec\n")
f.write(config['synthesis']['height']+"\n")
f.write("\n")
f.write("# Optical depth of the slab in the maximum of I (slab) or strength of the line (ME)\n")
if (nSlabs == 1):
	f.write(config['synthesis']['slab 1']['tau']+"\n")
else:
	if (includeff):
		f.write(config['synthesis']['slab 1']['tau']+"  "+config['synthesis']['slab 2']['tau']+"  "+config['synthesis']['ff']+"\n")
	else:
		f.write(config['synthesis']['slab 1']['tau']+"  "+config['synthesis']['slab 2']['tau']+"\n")		
f.write("\n")
f.write("# Source function enhancement\n") 
if (nSlabs == 1):
	f.write(config['synthesis']['slab 1']['beta']+"\n")
else:
	f.write(config['synthesis']['slab 1']['beta']+"  "+config['synthesis']['slab 2']['beta']+"\n")
f.write("\n")
f.write("# Boundary Stokes parameters (I0,Q0,U0,V0)  4.098093d-5 for 10830 A at disk center\n")
if (type(config['synthesis']['boundary condition']) is list):
	f.write("'single'\n")
	f.write("  ".join(config['synthesis']['boundary condition'])+"\n")
else:
	f.write("'file'\n")
	f.write("'{0}'\n".format(config['synthesis']['boundary condition']))
f.write("\n")
f.write("# Transition where to compute the emergent Stokes profiles\n")
f.write(multiplets[config['general parameters']['multiplet']]+"\n")
f.write("\n")
f.write("# Include atomic level polarization? (0-> no, 1-> yes)\n")
f.write(yesno[config['general parameters']['include atomic level polarization']]+"\n")
f.write("\n")
f.write("# Observation angle with respect to the local solar vertical theta,chi,gamma [degrees]\n")
f.write("  ".join(config['general parameters']['line-of-sight angles'])+"\n")
f.write("\n")
f.write("# Wavelength axis: minimum, maximum and number of grid points\n")
f.write("  ".join(config['general parameters']['wavelength axis'])+"\n")
f.write("\n")
f.write("# Line wavelength [A], Doppler velocity [km/s] and damping [a]\n")
if (nSlabs == 1):
	f.write(lambda0[config['general parameters']['multiplet']]+"  "+config['synthesis']['slab 1']['vdopp']+"  "+config['synthesis']['a']+"\n")
else:
	f.write(lambda0[config['general parameters']['multiplet']]+"  "+config['synthesis']['slab 1']['vdopp']+"  "+config['synthesis']['slab 2']['vdopp']+"  "+config['synthesis']['a']+"\n")
f.write("\n")
f.write("# Macroscopic velocity [km/s] (>0 is a redshift)\n")
if (nSlabs == 1):
	f.write(config['synthesis']['slab 1']['vmac']+"\n")
else:
	f.write(config['synthesis']['slab 1']['vmac']+"  "+config['synthesis']['slab 2']['vmac']+"\n")
f.write("\n")
f.write("# Include magneto-optical effects in the RT\n")
f.write(yesno[config['general parameters']['include magneto-optical effects in the rt']]+"\n")
f.write("\n")
f.write("# Include stimulated emission in the RT\n")
f.write(yesno[config['general parameters']['include stimulated emission in the rt']])
f.close()

# Run the code
args = ['mpiexec','-n',sys.argv[2],'./hazel','config_inversion.dat']

if (len(sys.argv) == 4):
	args.append(sys.argv[3])

if (len(sys.argv) == 5):
	args.append(sys.argv[3])
	args.append(sys.argv[4])

try:
	call(args)
except:
	print("A problem occured. Exiting...")
