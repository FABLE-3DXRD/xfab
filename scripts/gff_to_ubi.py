"""
Usage:
python gff_to_ubi input.gff output.ubi
"""

from ImageD11.columnfile import columnfile
from ImageD11.grain import grain, write_grain_file
import sys, numpy

c=columnfile(sys.argv[1])

labs = "UBI11 UBI12 UBI13 UBI21 UBI22 UBI23 UBI31 UBI32 UBI33".split()

grains = []

for i in range(c.nrows):
    ubi = [getattr(c, l)[i] for l in labs]
    t = c.x[i]*1000 , c.y[i]*1000, c.z[i]*1000
    grains.append( grain( numpy.reshape(ubi , (3,3)),
                          translation = t ) )


write_grain_file( sys.argv[2], grains )
