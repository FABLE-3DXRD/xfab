import numpy as n
import logging

from xfab import tools

def Uij2betaij(adp,ucell):
    """
    Uij2betaij transform the ADP U-matrix into the beta form 
    
    betaij = Uij2betaij(adp,unit_cell)
    
    INPUT:  adp: anisotropic displacement parameter U matrix
            Uijs given in order: [U11, U22, U33, U23, U13, U12]
            unit_cell = [a b c alpha beta gamma] 

    OUTPUT: betaij: beta displacement matrix
    
    Henning Osholm Sorensen, Risoe National Laboratory, June 23, 2006.
    Translated to python code March 29, 2008
    """
    
    U  = n.array([[adp[0], adp[5], adp[4]],
                  [adp[5], adp[1], adp[3]], 
                  [adp[4], adp[3], adp[2]]])

    betaij = n.zeros(3,3)
    cellstar = cellinvert(ucell)
 
    for i in range(3):
        for j in range(3):
            betaij[i,j] = 2*n.pi**2*cellstar[i]*cellstar[j]*U[i,j];


    return betaij

def StructureFactor(hkl,ucell,atomlist,sg,):
    """
     Calculation of the structure factor of reflection hkl
    
    [Freal Fimg] = StructureFactor(hkl,unit_cell,atomparam,sg,atom)
    
     INPUT : hkl =       [h, k, l] 
             unit_cell = [a, b, c, alpha, beta, gamma] 
             atomlist:  structural parameters
             sg:         read from sglib(sgno) (space group library)
             atom:       read from atomlib (library of form factors)
     OUTPUT: The real and imaginary parts of the the structure factor
    
     Henning Osholm Sorensen, June 23, 2006.
     Translated to python code March 29, 2008
     """
    
    stl = tools.sintl(ucell,hkl)
    noatoms = len(atomlist)

    Freal = 0.0
    Fimg = 0.0

    for i in range(noatoms):
        sfreal = 0.0
        sfimg = 0.0
    
        #Check whether isotrop or anisotropic displacements 
        Uelements = atomlist[i][5]

        if Uelements == 1:
            U = atomparam[i,6]
        elif  Uelements == 6:
            # transform Uij to betaij
            betaij = Uij2betaij(atomlist[i,6:12],cell);
        else:
            logging.error("wrong no of elements in atomlist")

        for j in range(sg.nosymop):
            # atomic displacement factor
            if Uelements == 1:
                expij=exp(-2*n.pi**2*U*stl**2);
            else:
                betaijrot = sg.rot[j]*betaij*sg.rot[j]
            # exponent for phase factor
            r = sg.rot[j]*atomlist[i,2:5] + sg.trans[j]
            exponent = 2*n.pi*hkl*r

            #forming the real and imaginary parts of F
            sfreal = sfreal + expij*n.cos(exponent)
            sfimg = sfimg + expij*n.sin(exponent)


        # Including the atomic formfactor
        formfac = FormFactor(atomlist[i,1],stl)*atomlist[12]/atomlist[13]
        Freal = Freal + formfac*sfreal
        Fimg = Fimg + formfac*sfimg

    return [Freal, Fimg]

def FormFactor(atom_no,stl):
    # Calculation of the atomic form factor at a specified sin(theta)/lambda
    # using the analytic fit to the Direc form factors from 
    # Int. Tab. Cryst Sect. C

    # Read atom library 
    f = open('atomlib.dat','r')
    data = f.readlines()[atom_no].split()
    f.close()

    # Calc form factor
    formfac = 0
    for i in range(4):
        formfac = formfac + eval(data[i+1])*n.exp(-eval(data[i+5])*stl*stl) 
    formfac = formfac + eval(data[9])
    return formfac

