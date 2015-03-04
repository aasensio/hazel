#!/usr/bin/env python
from configobj import ConfigObj
import sys
import os
from subprocess import call

if (len(sys.argv) < 2):
	print "Example usage: runHazel conf.ini"
	exit()

print "Using configuration file = "+sys.argv[1]

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

file = open('conf.input','w')