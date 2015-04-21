#!/usr/bin/env python


'''
tweakdetpars

This script optimizes the spatial detector-sample parameters so that
the number of peaks found is maximized. The peak-finding is applied
to a select number of seed grains (.ubi file), which is specified
as this script is called.

The script should be called like so:

tweakdetpars FF_fit.par m00.flt mall.ubi 0.02 2 FF_final.par -x FF
tweakdetpars NF_fit.par peaks_min10_t100.flt n01.ubi 0.03 3 NF_final.par -x NF

where 0.02 is the hkl tolerance and 2 is the number of iterations.

User input, such as parameter steps and range, can be found below the
function definitions.

'''

import argparse
import os
import shutil
import sys # where timeStep and layer are stored
import textwrap

import numpy as np
from scipy.optimize import leastsq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt # for plotting data if needed
import string

from ImageD11.refinegrains import parameters, refinegrains, transform


def scor_old(r): # find number of peaks
    '''
    Finds number of peaks for specified detector-sample parameter values.

    Parameters:
        r = parameter values
    Returns:
        n = number of peaks found
        e = ?
    '''

    #r.assignlabels() ## works here... try on line 202 instead

    npks = 0
    drlv = 0
    ks = r.grains.keys()
    for key in ks:
        g = r.grains[key]
        grainname = key[0]
        scanname = key[1]
        g.peaks_xyz = transform.compute_xyz_lab(
            [g.sc, g.fc],
            **r.parameterobj.parameters
            )
        # compute gv using current parameters, including grain position
        r.set_translation( key[0], key[1] )
        r.compute_gv(g)
        r.refine(g.ubi)
        npks += r.npks
        drlv += r.avg_drlv2
    return npks, drlv


def scor_maggie(r):
    # slightly modified version of refinegrains.gof
    drlv = 0.
    npks = 0.

    # defaulting to fitting all grains
    for key in r.grains.keys():

        g = r.grains[key]
        grainname = key[0]
        scanname = key[1]

        # Compute gv using current parameters
        # Keep labels fixed
        g.peaks_xyz = transform.compute_xyz_lab(
            [g.sc, g.fc],
            **r.parameterobj.parameters
            )
        r.compute_gv(g)

        g.set_ubi(r.refine(r.ubisread[grainname]))

        drlv += r.avg_drlv2
        npks += r.npks
    return npks, drlv


def scor_new(r):
    r.grains_to_refine = r.grains.keys()
    r.recompute_xlylzl = True
    return r.gof(), 0
    
def scor(r):
    """
    This scor function replaces all of the above since it should not refine grains (ubi and trans), but only assign reflections
    """
    r.assignlabels(quiet=True) 
    npks = np.sum(labels > -1 for labels in r.scandata[args.fltfile].labels)
    drlv = (np.sum(r.scandata[args.fltfile].drlv2)-(len(r.scandata[args.fltfile].drlv2)-npks))/npks
    return npks, drlv


def gaussian(x, *p):
    '''
    Creates y values for a specified Gaussian curve.

    Parameters:
        x = array of x values
        p* = array of Gaussian parameters (A, mu, sigma)
    Returns:
        y = array of y values for Gaussian curve
    '''
    A, mu, sigma = p
    y = A*np.exp(-(x-mu)**2/(2.*sigma**2))
    return y


def residuals(p, y, x): # difference between points and Gaussian curve
    A, mu, sigma = p
    err = y - A*np.exp(-(x-mu)**2/(2.*sigma**2))
    return err


def gauss_fit(x, y):
    '''
    Creates Gaussian fit of xy points.

    Parameters:
        x = array of x values (must be same # of elements as y!)
        y = array of y values (must be same # of elements as x!)
    Returns:
        p_opt = array of Gaussian parameters (A, mu, sigma)

    '''
    if y.std() == 0:
        # can't fit something that isn't varying
        return

    # set initial Gaussian parameters
    mu = x[np.argmax(y)] # set initial mu to x value of maximum
    A = np.amax(y) # set initial amplitude to y value of maximum
    # set initial sigma to x-distance from mu after 34.1% of curve area
    total_area = y.sum() # area under curve
    area_from_center = 0.
    sigma = 0
    for i in range( np.argmax(y), len(y) ):
        area_from_center += y[i]
        if area_from_center > .341 * total_area:
            sigma = x[i] - x[np.argmax(y)]
            break
    if sigma == 0:
        print('Range does not span a stddev, try increasing the range')
        sigma =  x[-1] - x[np.argmax(y)]
    p0 = [A, mu, sigma]
    print('Initial estimates of Gaussian fit parameters:')
    print('A     %.1f' % A)
    print('mu    %.4f' % mu)
    print('sigma %.4f' % sigma)

    # find optimal Gaussian parameters
    #p_opt = curve_fit(gaussian, x, y, p0=p0)
    return leastsq(residuals, p0, args=(y, x))[0]


def main(args):

    # set up peak-finding
    print
    r = refinegrains()
    r.loadparameters(args.parfile)
    r.loadfiltered(args.fltfile)
    r.readubis(args.ubifile)
    r.generate_grains()
    r.tolerance = float(args.tol)
    p = parameters.parameters()
    p.loadparameters(args.parfile)


    #*********************** BEGIN USER INPUT ***********************#

    print
    param_range = np.array([
    #     MIN          MAX
    #----------------------------
     [ r.parameterobj.parameters['distance'] * (1-args.dd),
       r.parameterobj.parameters['distance'] * (1+args.dd), ],     #dd
     [ r.parameterobj.parameters['y_center'] - args.dy,
       r.parameterobj.parameters['y_center'] + args.dy,     ],     #dy
     [ r.parameterobj.parameters['z_center'] - args.dz,
       r.parameterobj.parameters['z_center'] + args.dz,     ],     #dz
     [ r.parameterobj.parameters['tilt_x']   - args.dtx,
       r.parameterobj.parameters['tilt_x']   + args.dtx,    ],     #dtx
     [ r.parameterobj.parameters['tilt_y']   - args.dty,
       r.parameterobj.parameters['tilt_y']   + args.dtz,    ],     #dty
     [ r.parameterobj.parameters['tilt_z']   - args.dtz,
       r.parameterobj.parameters['tilt_z']   + args.dtz,    ],     #dtz
     [ r.parameterobj.parameters['wedge']    - args.dw,
       r.parameterobj.parameters['wedge']    + args.dw,     ]  ])  #dw

    #************************* END USER INPUT *************************#


    # read and record initial optimal paramter values
    init_pars = open(args.parfile,'r')
    init_pars_array = np.genfromtxt(init_pars,usecols=1)
    distance = init_pars_array[8]
    y_center = init_pars_array[25]
    z_center = init_pars_array[27]
    tilt_x = init_pars_array[20]
    tilt_y = init_pars_array[21]
    tilt_z = init_pars_array[22]
    wedge = init_pars_array[24]

    p_opt_array = [distance,y_center,z_center,tilt_x,tilt_y,tilt_z,wedge]

    # create list of parameter names

    all_param_names = ['distance', 'y_center', 'z_center',
                       'tilt_x', 'tilt_y', 'tilt_z', 'wedge']
    param_names = []
    if args.dd != 0:
        param_names.append('distance')
    if args.dy != 0:
        param_names.append('y_center')        
    if args.dz != 0:
        param_names.append('z_center')        
    if args.dtx != 0:
        param_names.append('tilt_x')        
    if args.dty != 0:
        param_names.append('tilt_y')        
    if args.dtz != 0:
        param_names.append('tilt_z')        
    if args.dw != 0 and (args.experiment == 'FF' or args.experiment == 'ff'):
        param_names.append('wedge')        


    print('varying:')
    print(param_names)
    print

    # find number of peaks using initial optimal parameter values
    for i, name in enumerate(all_param_names):
        r.parameterobj.parameters[name] = p_opt_array[i]
    r.assignlabels() ## try here...
    n, e = scor(r)
    print('Peaks found before parameter refinement: %d' % n)
    print

    # parameter refinement iteration loop
    for current_iteration in range(args.iterations):

        if current_iteration > 0:
        ## Update range on parameters
            param_range[0,0] = r.parameterobj.parameters['distance']*(1-args.dd)
            param_range[0,1] = r.parameterobj.parameters['distance']*(1+args.dd)
            param_range[1,0] = r.parameterobj.parameters['y_center']-args.dy
            param_range[1,1] = r.parameterobj.parameters['y_center']+args.dy
            param_range[2,0] = r.parameterobj.parameters['z_center']-args.dz
            param_range[2,1] = r.parameterobj.parameters['z_center']+args.dz
            param_range[3,0] = r.parameterobj.parameters['tilt_x']-args.dtx
            param_range[3,1] = r.parameterobj.parameters['tilt_x']+args.dtx
            param_range[4,0] = r.parameterobj.parameters['tilt_y']-args.dty
            param_range[4,1] = r.parameterobj.parameters['tilt_y']+args.dty
            param_range[5,0] = r.parameterobj.parameters['tilt_z']-args.dtz
            param_range[5,1] = r.parameterobj.parameters['tilt_z']+args.dtz
            param_range[6,0] = r.parameterobj.parameters['wedge']-args.dw
            param_range[6,1] = r.parameterobj.parameters['wedge']+args.dw
    # skip parameter refinement?
        if args.skip_par_ref:
            print('Option to skip parameter refinement is set to ON.'
                  'GrainSweeper will use parameters values from the nf_X.par'
                  'file, unless the opt_par_values_nf_X.txt already exists '
                  '(in which case values from this .txt file will be used.')
            print
            print('The parameter refinement skipping option can be '
                  'toggled/overridden in the tweakpars.py script.')
            print
            break

            ## set initial optimal parameters
            for i, name in enumerate(all_param_names):
                r.parameterobj.parameters[name] = p_opt_array[i]

        # create 3D array where:
        #     rows = parameter variation increments
        #     cols = parameter value, peaks found, e (3 cols total)
        #     slices = parameters varied
        param_peaks = np.zeros( shape= (args.steps, 3, len(param_names)) )

        ## find peaks for each par variation, fit Gaussian to curve, and find
        ## optimal par values
        for i, name in enumerate(param_names): # for each parameter
            print
            print('Varying parameter: %s; value: %g' % (name, p_opt_array[all_param_names.index(name)]))
            print('   distance   y_center  z_center    tilt_x      tilt_y'
                  '      tilt_z      wedge     peaks found     e')
            var_range = list(np.linspace(param_range[all_param_names.index(name),0], param_range[all_param_names.index(name),1],
                                         args.steps))
            for j, d in enumerate(var_range): # each parameter variation step
                r.parameterobj.parameters[name] = d
                dd = r.parameterobj.parameters['distance']
                yy = r.parameterobj.parameters['y_center']
                zz = r.parameterobj.parameters['z_center']
                tx = r.parameterobj.parameters['tilt_x']
                ty = r.parameterobj.parameters['tilt_y']
                tz = r.parameterobj.parameters['tilt_z']
                ww = r.parameterobj.parameters['wedge']

                n, e = scor(r)
                param_peaks[j,:,i] = d, n, e # write values to array
                print(
                    '%11.1f %9.1f %9.1f %11.6f %11.6f %11.6f %11.6f %10d '
                    '%11.5f' % (dd, yy, zz, tx, ty, tz, ww, n, e)
                    )

            # create Gaussian fit
            p_opt = gauss_fit( param_peaks[:,0,i], param_peaks[:,1,i] )
            if p_opt is None:
                print 'cant fit something that doesnt vary'
                r.parameterobj.parameters[name] = p_opt_array[all_param_names.index(name)]
                continue

            # create Gaussian xy values, and plot original and fit curves
            # together
            x_gauss_steps = np.linspace(param_range[all_param_names.index(name),0], param_range[all_param_names.index(name),1],
                                        200 )
            y_gauss_steps = gaussian(x_gauss_steps, *p_opt)
            plt.xlabel(str(name))
            plt.ylabel('peaks found')
            plt.plot(param_peaks[:,0,i], param_peaks[:,1,i])
            plt.plot(x_gauss_steps, y_gauss_steps)
            if args.experiment:
                figurename = ('%s_param_var_%s_iter%d.pdf'
                              % (args.experiment, name,
                                 current_iteration+1)
                              )
            else:
                figurename = ('NEW_param_var_%s_iter%d.pdf'
                              % (name, current_iteration+1)
                              )
                plt.savefig(figurename)
                plt.clf()

            ## print and record newest updated optimal parameter value,
            ## print peaks found
            print('Optimal parameter value: %.4f' % p_opt[1])
            p_opt_array[all_param_names.index(name)] = p_opt[1]
            r.parameterobj.parameters[name] = p_opt_array[all_param_names.index(name)]
            print
            n, e = scor(r) # find number of peaks using all optimal
                           #parameter values
            print('Peaks found after variation: %d' % n)

        # print optimal parameters from current iteration
        print
        print('Optimal parameter values for parameter refinement'
              'iteration %d:' % (current_iteration+1))
        print(' '.join(map(str, all_param_names)))
        print('  '.join(map(str, p_opt_array)))
        print

    # create optimal parameter values file (for importing into shell)
    if args.experiment:
        filename = '%s_opt_par_values.ini' % (args.experiment)
        writeline = args.experiment + '=' + str(n)
    else:
        filename = 'NEW_opt_par_values.ini'
        writeline = 'NEW = ' + str(n)

    text_file = open(filename, 'w')
    for i, name in enumerate(all_param_names):
        text_file.write( str(name) + '=' + str(p_opt_array[i]) + '\n')

    text_file.write(writeline)
    text_file.close()

    # write new paramterfile
    r.parameterobj.parameters['t_x'] = 0
    r.parameterobj.parameters['t_y'] = 0
    r.parameterobj.parameters['t_z'] = 0
    r.parameterobj.saveparameters(args.outparfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='tweakdetpars',
        description='Tweaks detector parameters to optimize reflections per '
                     'grain.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent('''
            examples:

            $ tweakdetpars FF_fit.par m00.flt mall.ubi 0.02 2 FF_final.par -x FF
            $ tweakdetpars NF_fit.par n00.flt nall.ubi 0.03 3 NF_final.par -x NF
            ''')
        )
    parser.add_argument(
        'parfile', type=str, help='name and path to starting detector par file'
        )
    parser.add_argument(
        'fltfile', type=str, help='name and path to flt file'
        )
    parser.add_argument(
        'ubifile', type=str, help='name and path to ubi file'
        )
    parser.add_argument(
        'tol', type=float, help='hkl tolerance'
        )
    parser.add_argument(
        'iterations', type=int, help='the number of iterations'
        )
    parser.add_argument(
        'outparfile', type=str, help='name and path to output detector par file'
        )
    parser.add_argument(
        '-s', '--steps', type=int, default=25,
        help='number of steps in each optimization calculation'
        )
    parser.add_argument(
        '-x', '--experiment', type=str, choices=['nf', 'ff', 'NF', 'FF'],
        default='', help='near-field or far-field experiment'
        )
    parser.add_argument(
        '-v', '--verbose', action='store_true',
        help='report progress in terminal'
        )
    parser.add_argument(
        '-k', '--skip-par-ref', action='store_true',
        help='skip the refinement process, uses the initial optimal par values'
        )
    parser.add_argument(
        '-dd', '--dd', type=float, default=0.05,
        help='range for tweaking parameter distance'
        )
    parser.add_argument(
        '-dy', '--dy', type=float, default=10,
        help='range for tweaking parameter y-center'
        )
    parser.add_argument(
        '-dz', '--dz', type=float, default=10,
        help='range for tweaking parameter z-center'
        )
    parser.add_argument(
        '-dtx', '--dtx', type=float, default=0.05,
        help='range for tweaking parameter tilt_x'
        )
    parser.add_argument(
        '-dty', '--dty', type=float, default=0.1,
        help='range for tweaking parameter tilt_y'
        )
    parser.add_argument(
        '-dtz', '--dtz', type=float, default=0.1,
        help='range for tweaking parameter tilt_z'
        )
    parser.add_argument(
        '-dw', '--dw', type=float, default=1,
        help='range for tweaking parameter wedge (not refined for NF)'
        )
    args = parser.parse_args()
    main(args)
