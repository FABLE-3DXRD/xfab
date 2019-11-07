#!/usr/bin/env python



##########################################################################
# IMPORT PYTHON PACKAGES
from __future__ import absolute_import
from __future__ import print_function
import numpy as np
import math
from xfab import tools
from ImageD11 import parameters
import os
import argparse

import textwrap
from six.moves import range
##########################################################################
myparms = parameters.parameters()

def readparfile(filename):

    myparms.loadparameters(filename)
    omat = np.zeros((2,2))
    if 'cell__a' in myparms.get_parameters().keys():
        a0 = float(myparms.get('cell__a'))
    if 'cell__b' in myparms.get_parameters().keys():
        b0 = float(myparms.get('cell__b'))
    if 'cell__c' in myparms.get_parameters().keys():
        c0 = float(myparms.get('cell__c'))
    if 'cell_alpha' in myparms.get_parameters().keys():
        alpha0 = float(myparms.get('cell_alpha'))
    if 'cell_beta' in myparms.get_parameters().keys():
        beta0 = float(myparms.get('cell_beta'))
    if 'cell_gamma' in myparms.get_parameters().keys():
        gamma0 = float(myparms.get('cell_gamma'))
    if 'cell_lattice_[P,A,B,C,I,F,R]' in myparms.get_parameters().keys():
        lattice = str(myparms.get('cell_lattice_[P,A,B,C,I,F,R]'))
    if 'chi' in myparms.get_parameters().keys():
        chi = float(myparms.get('chi'))
    if 'distance' in myparms.get_parameters().keys():
        distance = float(myparms.get('distance'))
    if 'fit_tolerance' in myparms.get_parameters().keys():
        fit_tolerance = float(myparms.get( 'fit_tolerance'))
    if 'min_bin_prob' in myparms.get_parameters().keys():
        min_bin_prob = float(myparms.get('min_bin_prob'))
    if 'no_bins' in myparms.get_parameters().keys():
        no_bins = float(myparms.get('no_bins'))
    if 'o11' in myparms.get_parameters().keys():
        omat[0,0] = float(myparms.get('o11'))
    if 'o12' in myparms.get_parameters().keys():
        omat[0,1] = float(myparms.get('o12'))
    if 'o21' in myparms.get_parameters().keys():
        omat[1,0] = float(myparms.get('o21'))
    if 'o22' in myparms.get_parameters().keys():
        omat[1,1] = float(myparms.get('o22'))
    if 'omegasign' in myparms.get_parameters().keys():
        omegasign = float(myparms.get('omegasign'))
    if 't_x' in myparms.get_parameters().keys():
        t_x = float(myparms.get('t_x'))
    if 't_y' in myparms.get_parameters().keys():
        t_y = float(myparms.get('t_y'))
    if 't_z' in myparms.get_parameters().keys():
        t_z = float(myparms.get('t_z'))
    if 'tilt_x' in myparms.get_parameters().keys():
        tilt_x = float(myparms.get('tilt_x'))
    if 'tilt_y' in myparms.get_parameters().keys():
        tilt_y = float(myparms.get('tilt_y'))
    if 'tilt_z' in myparms.get_parameters().keys():
        tilt_z = float(myparms.get('tilt_z'))
    if 'wavelength' in myparms.get_parameters().keys():
        wavelength = float(myparms.get('wavelength'))
    if 'wedge' in myparms.get_parameters().keys():
        wedge = float(myparms.get('wedge'))
    if 'y_center' in myparms.get_parameters().keys():
        y_center = float(myparms.get('y_center'))
    if 'y_size' in myparms.get_parameters().keys():
        y_size = float(myparms.get('y_size'))
    if 'z_center' in myparms.get_parameters().keys():
        z_center = float(myparms.get('z_center'))
    if 'z_size' in myparms.get_parameters().keys():
        z_size = float(myparms.get('y_size'))

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
    
    print('\nCreated files %s and %s.' % (sp1FLTfile,sp2FLTfile))


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
