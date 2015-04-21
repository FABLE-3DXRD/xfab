#!/usr/bin/env python
"""
Usage:
python gff_to_ubi input.gff output.ubi (detector.par)
NB! Pars only needed if gff does not contain ubis
"""

from ImageD11.columnfile import columnfile
from ImageD11.grain import grain, write_grain_file
import sys, numpy
from ImageD11 import parameters as ip
from xfab import tools

if len(sys.argv) < 3:
    print "\n"
    print "########################################################"
    print "Usage:"
    print "gff_to_ubi.py input.gff output.ubi (detector.par)"
    print "NB! Pars only needed if gff does not contain ubis"
    print "########################################################"
    print "\n"

c=columnfile(sys.argv[1])

try:
    labs = "UBI11 UBI12 UBI13 UBI21 UBI22 UBI23 UBI31 UBI32 UBI33".split()
    p = ip.parameters()
    p.loadparameters(sys.argv[3])
    uc = [p.parameters['cell__a'],p.parameters['cell__b'],p.parameters['cell__c'],
          p.parameters['cell_alpha'],p.parameters['cell_beta'],p.parameters['cell_gamma']]
    print uc
except:
    pass
    


grains = []

for i in range(c.nrows):
    t = c.x[i]*1000 , c.y[i]*1000, c.z[i]*1000
    try:
        ubi = [getattr(c, l)[i] for l in labs]
        grains.append( grain( numpy.reshape(ubi , (3,3)),
                          translation = t ) )
    except:
        print c.grainno[i]
        U = numpy.array([[c.U11[i],c.U12[i],c.U13[i]],
                         [c.U21[i],c.U22[i],c.U23[i]],
                         [c.U31[i],c.U32[i],c.U33[i]]])
        eps = numpy.array([c.eps11[i],c.eps12[i],c.eps13[i],c.eps22[i],c.eps23[i],c.eps33[i]])
        B = tools.epsilon_to_b(eps,uc)/2/numpy.pi
        UB = numpy.dot(U,B)
        ubi = numpy.linalg.inv(UB)
        grains.append( grain(ubi, translation = t ) )


write_grain_file( sys.argv[2], grains )
