#!/usr/bin/env python



##########################################################################
# IMPORT PYTHON PACKAGES
import numpy as np
import math
#import matplotlib.pyplot as plt

#from ImageD11 import transformer
from xfab import tools

import os
import argparse

#import shutil
#import sys
import textwrap

##########################################################################
def readparfile(filename):

    f=open(filename,'r')

    omat = np.zeros((2,2))

    for i in range(29):
        line=f.readline()
        if line.split()[0] == 'cell__a':
            a0 = float(line.split()[-1])
        if line.split()[0] == 'cell__b':
                b0 = float(line.split()[-1])
        if line.split()[0] == 'cell__c':
                c0 = float(line.split()[-1])
        if line.split()[0] == 'cell_alpha':
                alpha0 = float(line.split()[-1])
        if line.split()[0] == 'cell_beta':
                beta0 = float(line.split()[-1])
        if line.split()[0] == 'cell_gamma':
                gamma0 = float(line.split()[-1])
        if line.split()[0] == 'cell_lattice_[P,A,B,C,I,F,R]':
                lattice = str(line.split()[-1])
        if line.split()[0] == 'chi':
                chi = float(line.split()[-1])
        if line.split()[0] == 'distance':
                distance = float(line.split()[-1])
        if line.split()[0] == 'fit_tolerance':
                fit_tolerance = float(line.split()[-1])
        if line.split()[0] == 'min_bin_prob':
                min_bin_prob = float(line.split()[-1])
        if line.split()[0] == 'no_bins':
                no_bins = float(line.split()[-1])
        if line.split()[0] == 'o11':
                omat[0,0] = float(line.split()[-1])
        if line.split()[0] == 'o12':
                omat[0,1] = float(line.split()[-1])
        if line.split()[0] == 'o21':
                omat[1,0] = float(line.split()[-1])
        if line.split()[0] == 'o22':
                omat[1,1] = float(line.split()[-1])
        if line.split()[0] == 'omegasign':
                omegasign = float(line.split()[-1])
        if line.split()[0] == 't_x':
                t_x = float(line.split()[-1])
        if line.split()[0] == 't_y':
                t_y = float(line.split()[-1])
        if line.split()[0] == 't_z':
                t_z = float(line.split()[-1])
        if line.split()[0] == 'tilt_x':
                tilt_x = float(line.split()[-1])
        if line.split()[0] == 'tilt_y':
                tilt_y = float(line.split()[-1])
        if line.split()[0] == 'tilt_z':
                tilt_z = float(line.split()[-1])
        if line.split()[0] == 'wavelength':
                wavelength = float(line.split()[-1])
        if line.split()[0] == 'wedge':
                wedge = float(line.split()[-1])
        if line.split()[0] == 'y_center':
                y_center = float(line.split()[-1])
        if line.split()[0] == 'y_size':
                y_size = float(line.split()[-1])
        if line.split()[0] == 'z_center':
                z_center = float(line.split()[-1])
        if line.split()[0] == 'z_size':
                z_size = float(line.split()[-1])

    f.close()

    unitcell=np.array([a0, b0, c0, alpha0, beta0, gamma0])
    t_array = np.array([t_x, t_y, t_z])
    tilts = np.array([tilt_x, tilt_y, tilt_z])
    centers = np.array([y_center, z_center])
    pixsize = np.array([y_size, z_size])
    

    return(unitcell,lattice,distance,omat,tilts,wavelength,wedge,centers,pixsize)
##########################################################################
def sortFLT(FLTfile,sp1FLTfile,sp2FLTfile,sp1_sinth,sp2_sinth,distance,centers,pixsize,wavelength):

## # filename = FF_foundpeaks/00/peaks_t100.flt
## #  sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id
##   139.5461  650.8555  -180.0000  11  142.8182  139.5461  650.8555  1.4521  1.2420  0.0921  1.0000  0.0000  0.0000  1571.0000  234117.0000  187.0000  139  651  -180.0000  138  141  650  652  -180.0000  -180.0000  -650.8555  139.5461  1  0  3
##   217.2704  588.8578  -180.0000  17  210.1176  217.2704  588.8578  1.7047  1.2876  0.1872  1.0000  0.0000  -0.0000  3572.0000  874808.0000  389.0000  217  589  -180.0000  215  220  587  590  -180.0000  -180.0000  -588.8578  217.2704  1  0  5

    all_sinth = np.concatenate([sp1_sinth,sp2_sinth])
    sp1 = len(sp1_sinth)

    f=open(FLTfile,'r')
    fsp1=open(sp1FLTfile,'w')
    fsp2=open(sp2FLTfile,'w')


    for i in range(2):
        line=f.readline()
        fsp1.write(line)
        fsp2.write(line)
    titles = line.split()
    
    eof = 0
    while not(eof):
        line = f.readline()
        s = line.split()
        if len(s)>0:
       
            xdist = (float(s[1])-centers[0])*pixsize[0]
            ydist = (float(s[0])-centers[1])*pixsize[1]
            raddist = math.sqrt(xdist**2 + ydist**2)
            
            twth = math.atan(raddist/distance)
            sintwth = math.sin(twth/2)#/wavelength
            #print s[0],s[1],'\t',math.degrees(twth)
            
            error = abs(all_sinth-sintwth)
            minerror = min(error)

            idx = np.flatnonzero(error==minerror)
            if idx < sp1:
                #print 'SPECIES 1: Adding line to %s' % sp1FLTfile
                fsp1.write('%s' % line)
            else:
                #print 'SPECIES 2: Adding line to %s' % sp2FLTfile
                fsp2.write('%s' % line)
            #print sintwth,all_sinth[idx]

            ## save entire line to new file
        else:
            eof = 1

    f.close()
    fsp1.close()
    fsp2.close()
    
    print '\nCreated files %s and %s.' % (sp1FLTfile,sp2FLTfile)


    return()
##########################################################################
def initializeHKL(PARfile,sg):


    # Read files for information from both .log and .gve
    (unitcell,lattice,distance,omat,tilts,
              wavelength,wedge,centers,pixsize) = readparfile(PARfile)

    pixx = 2048
    pixy = pixx


    ## Calculate largest possible ring portion on detector
    ## units in mm from edge of detector
    xdist_a = (pixx - centers[0]) * pixsize[0]
    xdist_b = (centers[0]) * pixsize[0]
    ydist_a = (pixy - centers[1]) * pixsize[1]
    ydist_b = (centers[1]) * pixsize[1]
    
    ring_rad = np.zeros(( 4 ))
    ring_rad[0] = math.sqrt(xdist_a**2 + ydist_a**2)
    ring_rad[1] = math.sqrt(xdist_a**2 + ydist_b**2)
    ring_rad[2] = math.sqrt(xdist_b**2 + ydist_a**2)
    ring_rad[3] = math.sqrt(xdist_b**2 + ydist_b**2)

    max_rad = np.max(ring_rad)
    max2th = math.atan(max_rad/distance)

    ## Calculates min and max sin2th on detector
    sintlmin = 0
    sintlmax = np.sin(( max2th )/2)/wavelength
    sintlmax = sintlmax * 1.05 # over estimates a bit

    # Lattice information for unit cell
    hkls = tools.genhkl_all(unitcell,sintlmin,sintlmax,sgno=int(sg),output_stl=True)

    sinth = np.zeros((hkls.shape[0]))
    for i in range(0,hkls.shape[0]):
        sinth[i] = hkls[i,3]*wavelength
    sinth_unique = np.unique(sinth)
        
    return(sinth_unique)
##########################################################################
def main(args):

    if not os.path.exists(args.sp1parfile):
    	raise IOError('File %s does not exist.' % args.sp1parfile)
    if not os.path.exists(args.sp2parfile):
    	raise IOError('File %s does not exist.' % args.sp2parfile)
    if not os.path.exists(args.fltfile):
    	raise IOError('File %s does not exist.' % args.fltfile)

    ##########################################################################
    ## Calculate ring position for each species
    ## ...for species one
    sp1_sinth = initializeHKL(args.sp1parfile,args.sp1sg)
    ## ...for species one
    sp2_sinth = initializeHKL(args.sp2parfile,args.sp2sg)

    
    ##########################################################################
    ## Read FLT file and sort into two species
    sp1fltfile =  '%s_%s' % (args.sp1parfile.split('.')[0],args.fltfile)
    sp2fltfile =  '%s_%s' % (args.sp2parfile.split('.')[0],args.fltfile)

    if os.path.exists(sp1fltfile) and not args.force:
    	raise IOError('File %s already exists.' % sp1fltfile)
    if os.path.exists(sp2fltfile) and not args.force:
    	raise IOError('File %s already exists.' % sp2fltfile)


    (unitcell,lattice,distance,omat,tilts,
              wavelength,wedge,centers,pixsize) = readparfile(args.sp1parfile)
    sortFLT(args.fltfile,sp1fltfile,sp2fltfile,sp1_sinth,sp2_sinth,distance,centers,pixsize,wavelength)

    return()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='flt_split_phases.py',
        description='Sorts flt file for 2 species diffraction data.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
             examples:
 
             $ 
             $ flt_sort.py Cu_FF.par 225 W_FF.par 229 m00.flt
             ''')
         )
    parser.add_argument(
        'sp1parfile', type=str, help='name for first detpar file'
        )
    parser.add_argument(
        'sp1sg', type=int, help='space group corresponding to first file'
        )        
    parser.add_argument(
        'sp2parfile', type=str, help='name for second detpar file'
        )
    parser.add_argument(
        'sp2sg', type=int, help='space group corresponding to second file'
        )  
    parser.add_argument(
        'fltfile', type=str, help='flt file'
        )
    parser.add_argument(
        '-f', '--force', default=False, action='store_true',
        help='force removal of report file'
        )  
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report progress in terminal'
        )
    args = parser.parse_args()
    main(args)
