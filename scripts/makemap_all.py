#!/usr/bin/env python

"""
 Runs makemap.py multiple times, removes grains without enough peaks
 Should enter nine variables: nX.ubi, n(X-1).flt, nX.flt, initial cuts, final cuts,
 	stop limit, min number of peaks, full path to .par file, crystal symmetry
"""

from __future__ import absolute_import
from __future__ import print_function
import sys
import os
from string import split,digits,letters#,replace
from six.moves import range


# read and check input file
if len(sys.argv) != 10:
	print("Should enter nine variables:") 
	print("nX.ubi, n(X-1).flt, nX.flt, ")
	print("initial cuts, final cuts, stop limit, ")
	print("min number of peaks, full path to .par file, crystal symmetry \n") 
	print("Example:")
	print("makemap_all.py nall.ubi n00.flt rest.flt 0.03 0.02 5 100 NF_final.par cubic\n")
	sys.exit()
else:
	UBI = sys.argv[1]			# Input grains in .ubi format, overwritten by the program
	FLT_IN = sys.argv[2]		# Input peaks in .flt format
	FLT_OUT = sys.argv[3]		# Unassignted peaks in .flt format
	ICUTS=eval(sys.argv[4])		# Initial cuts, deviation from integer hkl, typically 0.02
	FCUTS=eval(sys.argv[5])		# Final cuts, deviation from integer hkl, typically 0.01
	LIMIT=eval(sys.argv[6])     # Minimum difference in remaining peaks after running makemap.py
	MINPEAKS=eval(sys.argv[7])  # Minimum number of peaks a seed grain needs to be retained in the final list
	PARS=sys.argv[8]            # Full path to the parameters file
	CRYSYMM=sys.argv[9]         # Crystal symmetry, e.g. "cubic"


# Run makemap loop 1 with first threshold cut ICUTS
print("Beginnning makemap loop 1 with cuts (-t %f)" %ICUTS)
command = "makemap.py -p %s -u %s -U %s -f %s -F %s -t %f -s %s > makemapTempOut.txt" %(PARS,UBI,UBI,FLT_IN,FLT_OUT,ICUTS,CRYSYMM)
os.system(command)

f = open("makemapTempOut.txt",'r')
lines = f.readlines()
f.close()

pks0 = eval(lines[-2].translate(None,letters))
pks1 = eval(lines[-1].translate(None,letters))
print("Initial number of peaks %i" %pks0)
print("Number of remaining peaks %i" %pks1)
	
new_assigned = LIMIT + 1
counter = 1
while new_assigned > LIMIT:	
	counter = counter + 1
	command = "makemap.py -p %s -u %s -U %s -f %s -F %s -t %f -s %s > makemapTempOut.txt" %(PARS,UBI,UBI,FLT_IN,FLT_OUT,ICUTS,CRYSYMM)
	os.system(command)
	f = open("makemapTempOut.txt",'r')
	lines = f.readlines()
	f.close()
	pks2 = eval(lines[-1].translate(None,letters))
	new_assigned = pks1 - pks2
	pks1 = pks2

print("makemap.py was run %i times" %(counter))


# Remove grains with less than MINPEAKS peaks
f = open(UBI,'r')
lines = f.readlines()
f.close()

breaklinenum = -1
for i in range(len(lines)):
	if "npks" in lines[i]:
		npks = eval(lines[i][1:].translate(None,letters))
		if npks < MINPEAKS:
			breaklinenum = i - 3
			break
if breaklinenum > 0:
	print("Removing grains with less than %i peaks" %MINPEAKS)
	print("Deleting from line %i on!" %breaklinenum)
	f = open(UBI,'w')
	for i in range(breaklinenum):
		f.write("%s" %lines[i])
	f.close()
elif breaklinenum == 0:
	print("All seed grains are bad.")
	print("Deleting %s and %s files!" %(UBI,FLT_OUT))
	os.remove(UBI)
	os.remove(FLT_OUT)
else:
	print("All grains have at least %i peaks" %MINPEAKS)		
			
# Run makemap loop 2 with first threshold cut FCUTS
if os.path.exists(UBI):
	print("Beginnning makemap loop 2 with cuts (-t %f)" %FCUTS)
	command = "makemap.py -p %s -u %s -U %s -f %s -F %s -t %f -s %s > makemapTempOut.txt" %(PARS,UBI,UBI,FLT_IN,FLT_OUT,FCUTS,CRYSYMM)
	os.system(command)

	f = open("makemapTempOut.txt",'r')
	lines = f.readlines()
	f.close()

	pks0 = eval(lines[-2].translate(None,letters))
	pks1 = eval(lines[-1].translate(None,letters))
	print("Initial number of peaks %i" %pks0)
	print("Number of remaining peaks %i" %pks1)
	
	new_assigned = LIMIT + 1
	counter = 1
	while new_assigned > LIMIT:	
		counter = counter + 1
		command = "makemap.py -p %s -u %s -U %s -f %s -F %s -t %f -s %s > makemapTempOut.txt" %(PARS,UBI,UBI,FLT_IN,FLT_OUT,FCUTS,CRYSYMM)
		os.system(command)
		f = open("makemapTempOut.txt",'r')
		lines = f.readlines()
		f.close()
		pks2 = eval(lines[-1].translate(None,letters))
		new_assigned = pks1 - pks2
		pks1 = pks2

	print("makemap.py was run %i more times" %(counter))

os.remove("makemapTempOut.txt")