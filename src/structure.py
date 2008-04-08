import numpy as n
import logging

from xfab import tools
from xfab import sg
from xfab import atomlib



def StructureFactor(hkl,ucell,sgname, atoms):
    """
    Calculation of the structure factor of reflection hkl
    
    [Freal Fimg] = StructureFactor(hkl,unit_cell,sg,atoms)
    
    INPUT : hkl =       [h, k, l] 
            unit_cell = [a, b, c, alpha, beta, gamma] 
            atomlist:  structural parameters
            sg:         read from sglib(sgno) (space group library)
            atom:       read from atomlib (library of form factors)
    OUTPUT: The real and imaginary parts of the the structure factor
    
    Henning Osholm Sorensen, June 23, 2006.
    Translated to python code March 29, 2008
    """
    
    mysg = sg.sg(sgname=sgname) 
    stl = tools.sintl(ucell,hkl)
    noatoms = len(atoms)

    Freal = 0.0
    Fimg = 0.0

    for i in range(noatoms):
        sfreal = 0.0
        sfimg = 0.0
    
        #Check whether isotrop or anisotropic displacements 
        if atoms[i].adp_type == 'Uiso':
            U = atoms[i].adp
            expij=n.exp(-8*n.pi**2*U*stl**2)
        elif atoms[i].adp_type == 'Uani':
            # transform Uij to betaij
            betaij = Uij2betaij(atoms[i].adp,ucell);
        else:
            logging.error("wrong no of elements in atomlist")

        for j in range(mysg.nsymop):
            # atomic displacement factor
            if atoms[i].adp_type == 'Uani':
                betaijrot = n.dot(mysg.rot[j],n.dot(betaij,mysg.rot[j]))
                expij=n.exp(-n.dot(hkl,n.dot(betaijrot,hkl)))
                
            # exponent for phase factor
            r = n.dot(mysg.rot[j],atoms[i].pos) + mysg.trans[j]
            exponent = 2*n.pi*n.dot(hkl,r)

            #forming the real and imaginary parts of F
            sfreal = sfreal + expij*n.cos(exponent)
            sfimg = sfimg + expij*n.sin(exponent)


        # Including the atomic formfactor
        formfac = FormFactor(atoms[i].atomtype,stl)*atoms[i].occ/atoms[i].symmulti
        
        Freal = Freal + formfac*sfreal
        Fimg = Fimg + formfac*sfimg

    return [Freal, Fimg]

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
    print adp
    U  = n.array([[adp[0], adp[5], adp[4]],
                  [adp[5], adp[1], adp[3]], 
                  [adp[4], adp[3], adp[2]]])

    betaij = n.zeros((3,3))
    cellstar = tools.CellInvert(ucell)
 
    for i in range(3):
        for j in range(3):
            betaij[i,j] = 2*n.pi**2*cellstar[i]*cellstar[j]*U[i,j];


    return betaij



def FormFactor(atomtype,stl):
    """
     Calculation of the atomic form factor at a specified sin(theta)/lambda
     using the analytic fit to the Direc form factors from 
     Int. Tab. Cryst Sect. C
    """

    # Read atom library 
    #f = open('atomlib.dat','r')
    #data = f.readlines()[atom_no].split()
    #f.close()
    print atomtype
    data = atomlib.Direc_ff[atomtype]
    print data
    # Calc form factor
    formfac = 0
    for i in range(4):
        formfac = formfac + data[i]*n.exp(-data[i+4]*stl*stl) 
    formfac = formfac + data[8]
    return formfac


class atom_entry:
    def __init__(self,label=None, atomtype=None, pos=None, adp_type=None, adp=None, occ=None, symmulti=None):
        self.label = label
        self.atomtype = atomtype
        self.pos = pos
        self.adp_type = adp_type
        self.adp = adp
        self.occ = occ
        self.symmulti = symmulti

class atomlist:
    def __init__(self, sgname=None, sgno=None, cell=None):
        self.sgname = sgname
        self.sgno = sgno
        self.cell = cell
        self.atom = []
    def add_atom(self,label=None, atomtype=None, pos=None, adp_type=None, adp=None, occ=None, symmulti=None):
        self.atom.append(atom_entry(label=label, atomtype=atomtype, pos=pos, adp_type=adp_type,
                                    adp=adp, occ=occ, symmulti=symmulti))

class build_atomlist:
    def __init__(self):
        self.atomlist = atomlist()
        
    def CIFopen(self,ciffile=None,cifblkname=None):
        from CifFile import ReadCif # part of the PycifRW module
        try:
            cf = ReadCif(ciffile)
        except:
            print 'File %s could not be accessed' %ciffile

        if cifblkname == None:   
            #Try to guess blockname                                                     
            blocks = cf.keys()
            print blocks
            if len(blocks) > 1:
                if len(blocks) == 2 and 'global' in blocks:
                    blockname = blocks[abs(blocks.index('global')-1)]
                    print blockname
                else:
                    print 'More than one possible data set:'
                    print 'The following data block names are in the file:'
                    for block in blocks:
                        print block
                    raise IOError
            else:
                blockname = blocks[0]
        try:
            self.cifblk = cf[cifblkname]
        except:
            print 'Block - %s - not found in %s' %(blockname,ciffile)
            raise IOError
        return self.cifblk

    def CIFread(self,ciffile = None,cifblk = None, cifblkname = None):
        from re import sub
        from string import upper
        if cifblk == None:
            print 'in here'
            try:
                cifblk = self.CIFopen(ciffile=ciffile,cifblkname=cifblkname)
            except:
                raise IOError()
        print cifblk['_cell_length_a']
        self.atomlist.cell = [self.remove_esd(cifblk['_cell_length_a']),
                              self.remove_esd(cifblk['_cell_length_b']),
                              self.remove_esd(cifblk['_cell_length_c']),
                              self.remove_esd(cifblk['_cell_angle_alpha']),
                              self.remove_esd(cifblk['_cell_angle_beta']),
                              self.remove_esd(cifblk['_cell_angle_gamma'])]

        #self.atomlist.sgname = upper(sub("\s+","",cifblk['_symmetry_space_group_name_H-M']))
        self.atomlist.sgname = sub("\s+","",cifblk['_symmetry_space_group_name_H-M'])

        for i in range(len(cifblk['_atom_site_type_symbol'])):
            label = cifblk['_atom_site_label'][i]
            #atomno = atomtype[upper(cifblk['_atom_site_type_symbol'][i])]
            atomtype = upper(cifblk['_atom_site_type_symbol'][i])
            print atomtype
            x = self.remove_esd(cifblk['_atom_site_fract_x'][i])
            y = self.remove_esd(cifblk['_atom_site_fract_y'][i])
            z = self.remove_esd(cifblk['_atom_site_fract_z'][i])
            try:
                adp_type = cifblk['_atom_site_adp_type'][i]
            except:
                adp_type = None
            try:
                occ = self.remove_esd(cifblk['_atom_site_occupancy'][i])
            except:
                occ = 1.0

            if cifblk.has_key('_atom_site_symmetry_multiplicity'):
                multi = self.remove_esd(cifblk['_atom_site_symmetry_multiplicity'][i])
            # In old SHELXL versions this code was written as '_atom_site_symetry_multiplicity'
            elif cifblk.has_key('_atom_site_symetry_multiplicity'):
                multi = self.remove_esd(cifblk['_atom_site_symetry_multiplicity'][i])
            else:
                multi = 1.0
                print 'WARNING %s has no multiplicity given in cif file >>>  Value default 1.0' %label



            if adp_type == None:
                adp = 0.0
            elif adp_type == 'Uiso':
                adp = self.remove_esd(cifblk['_atom_site_U_iso_or_equiv'][i])
            elif adp_type == 'Uani':
                anisonumber = cifblk['_atom_site_aniso_label'].index(label)
                adp = [ self.remove_esd(cifblk['_atom_site_aniso_U_11'][anisonumber]),
                        self.remove_esd(cifblk['_atom_site_aniso_U_22'][anisonumber]),
                        self.remove_esd(cifblk['_atom_site_aniso_U_33'][anisonumber]),
                        self.remove_esd(cifblk['_atom_site_aniso_U_23'][anisonumber]),
                        self.remove_esd(cifblk['_atom_site_aniso_U_13'][anisonumber]),
                        self.remove_esd(cifblk['_atom_site_aniso_U_12'][anisonumber])]
            
            self.atomlist.add_atom(label=label, atomtype=atomtype, pos = [x,y,z],
                                   adp_type= adp_type, adp = adp, occ=occ , symmulti=multi)

    def remove_esd(self,a):
        """                                                                         
        This function will remove the esd part of the entry,
        fx '1.234(56)' to '1.234'.                                                                             
        """
        from string import atof
        
        if a.find('(') == -1:
            value = atof(a)
        else:
            value = atof(a[:a.find('(')])
        return value


atomtype = {'H'  :  1, 'HE' :  2, 'LI' :  3, 'BE' :  4, 'B'  :  5, 'C'  :  6,
            'N'  :  7, 'O'  :  8, 'F'  :  9, 'NE' : 10, 'NA' : 11, 'MG' : 12,
            'AL' : 13, 'SI' : 14, 'P'  : 15, 'S'  : 16, 'CL' : 17, 'AR' : 18,
            'K'  : 19, 'CA' : 20, 'SC' : 21, 'TI' : 22, 'V'  : 23, 'CR' : 24,
            'MN' : 25, 'FE' : 26, 'CO' : 27, 'NI' : 28, 'CU' : 29, 'ZN' : 30,
            'GA' : 31, 'GE' : 32, 'AS' : 33, 'SE' : 34, 'BR' : 35, 'KR' : 36,
            'RB' : 37, 'SR' : 38, 'Y'  : 39, 'ZR' : 40, 'NB' : 41, 'MO' : 42,
            'TC' : 43, 'RU' : 44, 'RH' : 45, 'PD' : 46, 'AG' : 47, 'CD' : 48,
            'IN' : 49, 'SN' : 50, 'SB' : 51, 'TE' : 52, 'I'  : 53, 'XE' : 54,
            'CS' : 55, 'BA' : 56, 'LA' : 57, 'CE' : 58, 'PR' : 59, 'PM' : 60,
            'SM' : 61, 'EU' : 62, 'GD' : 63, 'GD' : 64, 'TB' : 65, 'DY' : 66,
            'HO' : 67, 'ER' : 68, 'TM' : 69, 'YB' : 70, 'LU' : 71, 'HF' : 72,
            'TA' : 73, 'W'  : 74, 'RE' : 75, 'OS' : 76, 'IR' : 77, 'PT' : 78,
            'AU' : 79, 'HG' : 80, 'TL' : 81, 'PB' : 82, 'BI' : 83, 'PO' : 84,
            'AT' : 85, 'RN' : 86, 'FR' : 87, 'RA' : 88, 'AC' : 89, 'TH' : 90,
            'PA' : 91, 'U'  : 92, 'NP' : 93, 'PU' : 94, 'AM' : 95, 'CM' : 96}

