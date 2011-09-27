"""
Usage:
python ubi_to_gff.py input.ubi detector.par output.gff
"""

from ImageD11 import columnfile as ic
from ImageD11 import grain as ig
from ImageD11 import parameters as ip
import sys
import numpy as n
from string import split
from xfab import tools

list_of_grains = ig.read_grain_file(sys.argv[1])
p = ip.parameters()
p.loadparameters(sys.argv[2])
uc = [p.parameters['cell__a'],p.parameters['cell__b'],p.parameters['cell__c'],p.parameters['cell_alpha'],p.parameters['cell_beta'],p.parameters['cell_gamma']]
print uc

grainno = []
x = []
y = []
z = []
rodx = []
rody = []
rodz = []
U11 = []
U12 = []
U13 = []
U21 = []
U22 = []
U23 = []
U31 = []
U32 = []
U33 = []
eps11 = []
eps22 = []
eps33 = []
eps23 = []
eps13 = []
eps12 = []
titles = ["grainno","x","y","z","rodx","rody","rodz","U11","U12","U13","U21","U22","U23","U31","U32","U33","eps11","eps22","eps33","eps23","eps13","eps12"]
for i in range(len(list_of_grains)):
    grainno.append(eval(split(list_of_grains[i].name,':')[0]))
    x.append(list_of_grains[i].translation[0]/1000.)
    y.append(list_of_grains[i].translation[1]/1000.)
    z.append(list_of_grains[i].translation[2]/1000.)
    ubi = list_of_grains[i].ubi
    (U,eps) = tools.ubi_to_u_and_eps(ubi,uc)
    rod = tools.u_to_rod(U)
    rodx.append(rod[0])
    rody.append(rod[1])
    rodz.append(rod[2])
    U11.append(U[0,0])
    U12.append(U[0,1])
    U13.append(U[0,2])
    U21.append(U[1,0])
    U22.append(U[1,1])
    U23.append(U[1,2])
    U31.append(U[2,0])
    U32.append(U[2,1])
    U33.append(U[2,2])
    eps11.append(eps[0])
    eps12.append(eps[1])
    eps13.append(eps[2])
    eps22.append(eps[3])
    eps23.append(eps[4])
    eps33.append(eps[5])

gff = ic.newcolumnfile(titles)
gff.ncols = len(titles)
gff.nrows = len(grainno)
gff.bigarray = n.zeros((gff.ncols,gff.nrows))
gff.set_attributes()
gff.addcolumn(grainno,"grainno")
gff.addcolumn(x,"x")
gff.addcolumn(y,"y")
gff.addcolumn(z,"z")
gff.addcolumn(rodx,"rodx")
gff.addcolumn(rody,"rody")
gff.addcolumn(rodz,"rodz")
gff.addcolumn(U11,"U11")
gff.addcolumn(U12,"U12")
gff.addcolumn(U13,"U13")
gff.addcolumn(U21,"U21")
gff.addcolumn(U22,"U22")
gff.addcolumn(U23,"U23")
gff.addcolumn(U31,"U31")
gff.addcolumn(U32,"U32")
gff.addcolumn(U33,"U33")
gff.addcolumn(eps11,"eps11")
gff.addcolumn(eps22,"eps22")
gff.addcolumn(eps33,"eps33")
gff.addcolumn(eps23,"eps23")
gff.addcolumn(eps13,"eps13")
gff.addcolumn(eps12,"eps12")
gff.writefile(sys.argv[3])
