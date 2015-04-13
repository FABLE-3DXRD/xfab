#!/usr/bin/python

"""
tweakpars_NF.py

This script optimizes the spatial detector-sample parameters so that
the number of peaks found is maximized. The peak-finding is applied
to a select number of seed grains (.ubi file), which is specified
as this script is called.

The script should be called like so:

tweakpars.py NF_input.par peaks.flt grains.ubi 0.03 3 NF_output.par

where 0.03 is the hkl tolerance and 3 is the number of iterations.

NB! Note that the wedge is not fitted for the nearfield detector.

User input, such as parameter steps and range, can be found below the
function definitions.

"""

# import modules
print
import numpy as np
from scipy.optimize import leastsq
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt # for plotting data if needed
from ImageD11.refinegrains import * # located in ~/../../usr/lib64/python2.6/site-packages/ImageD11/
#from opt_par_values import * # import initial optimal parameter values
import sys # where timeStep and layer are stored
import string


# define functions:

def scor(r): # find number of peaks
    """
    Finds number of peaks for specified detector-sample parameter values.

    Parameters:
        r = parameter values
    Returns:
        n = number of peaks found
        e = ?
    """
    npks = 0
    drlv = 0
    ks = r.grains.keys()
    for key in ks:
        g = r.grains[key]
        grainname = key[0]
        scanname = key[1]
        # compute gv using current parameters, including grain position
        r.set_translation( key[0], key[1] )
        r.compute_gv(g)
        r.refine(g.ubi)
        npks += r.npks
        drlv += r.avg_drlv2
    return npks, drlv


def gaussian(x, *p):
    """
    Creates y values for a specified Gaussian curve.

    Parameters:
        x = array of x values
        p* = array of Gaussian parameters (A, mu, sigma)
    Returns:
        y = array of y values for Gaussian curve
    """
    A, mu, sigma = p
    y = A*np.exp(-(x-mu)**2/(2.*sigma**2))
    return y


def gauss_fit(x, y):
    """
    Creates Gaussian fit of xy points. Requires residuals() function.

    Parameters:
        x = array of x values (must be same # of elements as y!)
        y = array of y values (must be same # of elements as x!)
    Returns:
        p_opt = array of Gaussian parameters (A, mu, sigma)

    """
    def residuals(p, y, x): # difference between points and Gaussian curve
        A, mu, sigma = p
        err = y - A*np.exp(-(x-mu)**2/(2.*sigma**2))
        return err
    # set initial Gaussian parameters
    mu_int = x[np.argmax(y)] # set initial mu to x value of maximum
    A_int = np.amax(y) # set initial amplitude to y value of maximum
    # set initial sigma to x-distance from mu after 34.1% of curve area
    total_area_bins = np.zeros( len(y) )
    for i in range( 0, len(y) ):
        total_area_bins[i] = x[i] * y[i] # area bins under curve
    area_bins_from_center = 0.
    for i in range( np.argmax(y), len(y) ):
        area_bins_from_center += total_area_bins[i]
        if area_bins_from_center > .341 * sum(total_area_bins):
            sigma_int = x[i] - x[np.argmax(y)]
            break
    # if sigma can't be estimated, set it to delta x
    try: sigma_int
    except: sigma_int = x[1] - x[0]
    if sigma_int == 0:
        sigma_int = x[1] - x[0]
    p_int = [A_int, mu_int, sigma_int]
    print 'Initial estimates of Gaussian fit parameters:'
    print 'A     ', '%.1f' % A_int
    print 'mu    ', '%.4f' % mu_int
    print 'sigma ', '%.4f' % sigma_int
    #sigma:, '%.1f' % A_int, ' %.1f' % mu_int, ' %.1f' % sigma_int
    # find optimal Gaussian parameters
    #p_opt = curve_fit(gaussian, x, y, p0=p_int)
    p_opt = leastsq( residuals, p_int, args=(y, x) )
    p_opt = p_opt[0]
    return p_opt


# unpack input arguments
parfile, fltfile, ubifile, tol, iterations, outparfile = sys.argv[1:]
skip_par_ref = 'N'
tol = float(tol)
iterations=int(iterations)

# set up peak-finding
print
r = refinegrains()
r.loadparameters(parfile)
r.loadfiltered(fltfile)
r.readubis(ubifile)
r.generate_grains()
r.tolerance = float(tol)
p = parameters.parameters()
p.loadparameters(parfile)




#********************************** BEGIN USER INPUT **********************************#

# set number of parameter variation increments
param_steps = 25
params_to_vary = 6

# force skip the refinement process and just use the initial optimal par values?
#skip_par_ref='Y'
#skip_par_ref='N'

# force number of refinement iterations?
#iterations = 2

print
param_range = np.array([
#     MIN          MAX
#----------------------------
 [ r.parameterobj.parameters['distance']*0.95,  r.parameterobj.parameters['distance']*1.05,     ],     #dd
 [ r.parameterobj.parameters['y_center']-10.,    r.parameterobj.parameters['y_center']+10.,     ],     #dy
 [ r.parameterobj.parameters['z_center']-10.,    r.parameterobj.parameters['z_center']+10.,     ],     #dz
 [ r.parameterobj.parameters['tilt_x']-0.05,      r.parameterobj.parameters['tilt_x']+0.05,     ],     #dtx
 [ r.parameterobj.parameters['tilt_y']-0.1,      r.parameterobj.parameters['tilt_y']+0.1,     ],     #dty
 [ r.parameterobj.parameters['tilt_z']-0.1,      r.parameterobj.parameters['tilt_z']+0.1,     ],     #dtz
 [ r.parameterobj.parameters['wedge']-1.,        r.parameterobj.parameters['wedge']+1.,     ]  ]) #dw

#********************************** END USER INPUT **********************************#


# read and record initial optimal paramter values
init_pars = open(parfile,'r')
init_pars_array = np.genfromtxt(init_pars,usecols=1)
distance = init_pars_array[8]
y_center = init_pars_array[25]
z_center = init_pars_array[27]
tilt_x = init_pars_array[20]
tilt_y = init_pars_array[21]
tilt_z = init_pars_array[22]
wedge = init_pars_array[24]

p_opt_array = [distance,y_center,z_center,tilt_x,tilt_y,tilt_z,wedge]
# refine in different order?
#p_opt_array = [tilt_x,tilt_y,tilt_z,wedge,y_center,z_center,distance]

# create list of parameter names
#param_names = ['distance', 'y_center', 'z_center', 'tilt_x', 'tilt_y', 'tilt_z', 'wedge']
all_param_names = ['distance','y_center','z_center','tilt_x','tilt_y','tilt_z','wedge']
param_names = all_param_names[0:params_to_vary]
# refine in different order?
#param_names = ['tilt_x','tilt_y','tilt_z','wedge','y_center','z_center','distance']
print
print 'varying:' 
print param_names
print

# find number of peaks using initial optimal parameter values
for i, name in enumerate(all_param_names):
    r.parameterobj.parameters[name] = p_opt_array[i]
print
n, e = scor(r)
print 'Peaks found before parameter refinement:', n
print

# parameter refinement iteration loop
for current_iteration in range(iterations):

    if current_iteration > 0:
        param_range[0,0] = r.parameterobj.parameters['distance']*0.975
        param_range[0,1] = r.parameterobj.parameters['distance']*1.025 #changee range on dd
        param_range[1,0] = r.parameterobj.parameters['y_center']-10.
        param_range[1,1] = r.parameterobj.parameters['y_center']+10. #update range on dy
        param_range[2,0] = r.parameterobj.parameters['z_center']-10.
        param_range[2,1] = r.parameterobj.parameters['z_center']+10. #update range on dz
        param_range[3,0] = r.parameterobj.parameters['tilt_x']-0.05
        param_range[3,1] = r.parameterobj.parameters['tilt_x']+0.05 #update range on dtx
        param_range[4,0] = r.parameterobj.parameters['tilt_y']-0.1
        param_range[4,1] = r.parameterobj.parameters['tilt_y']+0.1 #update range on dty
        param_range[5,0] = r.parameterobj.parameters['tilt_z']-0.1
        param_range[5,1] = r.parameterobj.parameters['tilt_z']+0.1 #update range on dtz
        param_range[6,0] = r.parameterobj.parameters['wedge']-1.
        param_range[6,1] = r.parameterobj.parameters['wedge']+1. #update range on dw
    # skip parameter refinement?
    if skip_par_ref == 'Y':
        print("Option to skip parameter refinement is set to ON. GrainSweeper will use "
              "parameters values from the nf_X.par file, unless the opt_par_values_nf_X.txt "
              "already exists (in which case values from this .txt file will be used.")
        print
        print("The parameter refinement skipping option can be toggled/overridden in the tweakpars.py script.")
        print
        break

    # set initial optimal parameters
    for i, name in enumerate(all_param_names):
        r.parameterobj.parameters[name] = p_opt_array[i]

    # create 3D array where:
    #     rows = parameter variation increments
    #     cols = parameter value, peaks found, e (3 cols total)
    #     slices = parameters varied
    param_peaks = np.zeros( shape= (param_steps, 3, len(param_names)) )

    # find peaks for each par variation, fit Gaussian to curve, and find optimal par values
    for i, name in enumerate(param_names): # for each parameter
        print
        print 'Varying parameter:', name
        print '   distance   y_center  z_center    tilt_x      tilt_y      tilt_z      wedge     peaks found     e'
        var_range = list( np.linspace( param_range[i,0], param_range[i,1], param_steps ) )
        for j, d in enumerate(var_range): # for each parameter variation step
            r.parameterobj.parameters[name] = d
            dd = r.parameterobj.parameters['distance']
            yy = r.parameterobj.parameters['y_center']
            zz = r.parameterobj.parameters['z_center']
            tx = r.parameterobj.parameters['tilt_x']
            ty = r.parameterobj.parameters['tilt_y']
            tz = r.parameterobj.parameters['tilt_z']
            ww = r.parameterobj.parameters['wedge']
            n, e = scor(r) # find number of peaks
            param_peaks[j,:,i] = d, n, e # write values to 4D array
            print '%11.1f' % dd, '%9.1f' % yy, '%9.1f' % zz, '%11.6f' % tx, '%11.6f' % ty, '%11.6f' % tz, '%11.6f' % ww, '%10d' % n, '%11.5f' % e

        # create Gaussian fit
        p_opt = gauss_fit( param_peaks[:,0,i], param_peaks[:,1,i] )

        # create Gaussian xy values, and plot original and fit curves together
        x_gauss_steps = np.linspace( param_range[i,0], param_range[i,1], 200 )
        y_gauss_steps = gaussian(x_gauss_steps, *p_opt)
        plt.xlabel(str(name))
        plt.ylabel('peaks found')
        plt.plot(param_peaks[:,0,i], param_peaks[:,1,i])
        plt.plot(x_gauss_steps, y_gauss_steps)
        plt.savefig('NF_param_var_%s_iter%d.pdf' % (name, current_iteration+1) )
        plt.clf()

        # print and record newest updated optimal parameter value, print peaks found
        print 'Optimal parameter value:', '%.4f' % p_opt[1]
        p_opt_array[i] = p_opt[1]
        r.parameterobj.parameters[name] = p_opt_array[i]
        print
        n, e = scor(r) # find number of peaks using all optimal parameter values
        print 'Peaks found after variation:', n

    # print optimal parameters from current iteration
    print
    print 'Optimal parameter values for parameter refinement iteration %d:' % (current_iteration+1)
    print ' '.join(map(str, param_names))
    print '  '.join(map(str, p_opt_array))
    print

# create optimal parameter values file (for importing into shell)
text_file = open("NF_opt_par_values.ini", "w")
for i, name in enumerate(all_param_names):
    text_file.write( str(name) + "=" + str(p_opt_array[i]) + "\n")

text_file.write( "NF" + "=" + str(n) )

text_file.close()

# write new paramterfile
r.parameterobj.saveparameters(outparfile)

