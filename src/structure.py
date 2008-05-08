import numpy as n
import logging

from xfab import tools
from xfab import sg
from xfab import atomlib

def StructureFactor(hkl,ucell,sgname,atoms,disper = None):
    """
    Calculation of the structure factor of reflection hkl
    
    [Freal Fimg] = StructureFactor(hkl,unit_cell,sg,atoms)
    
    INPUT : hkl =       [h, k, l] 
            unit_cell = [a, b, c, alpha, beta, gamma] 
            sgname:     space group name (e.g. 'P 21/c')
            atoms:      structural parameters (as an object)
    OUTPUT: The real and imaginary parts of the the structure factor
    
    Henning Osholm Sorensen, June 23, 2006.
    Translated to python code April 8, 2008
    """
    
    mysg = sg.sg(sgname=sgname) 
    stl = tools.sintl(ucell,hkl)
    noatoms = len(atoms)

    Freal = 0.0
    Fimg = 0.0

    for i in range(noatoms):
        #Check whether isotrop or anisotropic displacements 
        if atoms[i].adp_type == 'Uiso':
            U = atoms[i].adp
            expij=n.exp(-8*n.pi**2*U*stl**2)
        elif atoms[i].adp_type == 'Uani':
            # transform Uij to betaij
            betaij = Uij2betaij(atoms[i].adp,ucell);
        else:
            logging.error("wrong no of elements in atomlist")

        # Atomic form factors
        f = FormFactor(atoms[i].atomtype,stl)
        if disper != None:
            fp = disper[atoms[i].atomtype][0]
            fpp = disper[atoms[i].atomtype][1]
        else:
            fp = 0.0
            fpp = 0.0

        for j in range(mysg.nsymop):
            # atomic displacement factor
            if atoms[i].adp_type == 'Uani':
                betaijrot = n.dot(mysg.rot[j],n.dot(betaij,mysg.rot[j]))
                expij=n.exp(-n.dot(hkl,n.dot(betaijrot,hkl)))
                
            # exponent for phase factor
            r = n.dot(mysg.rot[j],atoms[i].pos) + mysg.trans[j]
            exponent = 2*n.pi*n.dot(hkl,r)

            #forming the real and imaginary parts of F
            s = n.sin(exponent)
            c = n.cos(exponent)
            site_pop = atoms[i].occ/atoms[i].symmulti
            Freal = Freal + expij*(c*(f+fp)-s*fpp)*site_pop
            Fimg = Fimg + expij*(s*(f+fp)+c*fpp)*site_pop

    return [Freal, Fimg]

def Uij2betaij(adp,ucell):
    """
    Uij2betaij transform the ADP U-matrix into the beta form 
    
    betaij = Uij2betaij(adp,unit_cell)
    
    INPUT:  adp: anisotropic displacement parameter U matrix
            Uijs are given in this order: [U11, U22, U33, U23, U13, U12]
            unit_cell = [a, b, c, alpha, beta, gamma] 

    OUTPUT: betaij: beta displacement matrix
    
    Henning Osholm Sorensen, Risoe National Laboratory, June 23, 2006.
    Translated to python code March 29, 2008
    """
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
     using the analytic fit to the  form factors from 
     Int. Tab. Cryst Sect. C 6.1.1.4

     INPUT:   atomtype: Atom type (string) e.g. 'C' 
              stl: form factor calculated sin(theta)/lambda = stl
     OUTPUT:  atomic form factor (no dispersion)
              

     Henning Osholm Sorensen, Risoe National Laboratory, April 9, 2008.

    """

    data = atomlib.formfactor[atomtype]

    # Calc form factor
    formfac = 0

    for i in range(4):
        formfac = formfac + data[i]*n.exp(-data[i+4]*stl*stl) 
    formfac = formfac + data[8]

    return formfac


def int_intensity(F2,L,P,I0,wavelength,cell_vol,cryst_vol):
	"""
	Calculate the reflection intensities scaling factor
        
        INPUT:
        F2        : the structure factor squared
        L         : Lorentz factor
        P         : Polarisation factor
        I0        : Incoming beam flux
        wavelength: in Angstroem
        cell_vol  : Volume of unit cell in AA^3
        cryst_vol : Volume of crystal in mm^3

        OUTPUT:
        int_intensity: integrated intensity

        """
#        print F2,L,P,I0,wavelength,cell_vol,cryst_vol
        
	emass =9.1093826e-31
	echarge = 1.60217653e-19
	pi4eps0 = 1.11265e-10
	c = 299792458.0
	k1 = (echarge**2/(pi4eps0*emass*c**2)*1000)**2 # Unit is mm
	k2 = wavelength**3 * cryst_vol * 1e21/cell_vol**2 # 1e21 to go from mm^3 to AA^3
        return k1*k2*I0*L*P*F2


class atom_entry:
    def __init__(self,label=None, atomtype=None, pos=None,
                 adp_type=None, adp=None, occ=None, symmulti=None):
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
        self.dispersion = {}
        self.atom = []
    def add_atom(self,label=None, atomtype=None, pos=None, 
                 adp_type=None, adp=None, occ=None, symmulti=None):
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
            logging.error('File %s could not be accessed' %ciffile)

        if cifblkname == None:   
            #Try to guess blockname                                                     
            blocks = cf.keys()
            if len(blocks) > 1:
                if len(blocks) == 2 and 'global' in blocks:
                    cifblkname = blocks[abs(blocks.index('global')-1)]
                else:
                    logging.error('More than one possible data set:')
                    logging.error('The following data block names are in the file:')
                    for block in blocks:
                        logging.error(block)
                    raise Exception
            else:
                # Only one available
                cifblkname = blocks[0]
        #Extract block
        try:
            self.cifblk = cf[cifblkname]
        except:
            logging.error('Block - %s - not found in %s' %(blockname,ciffile))
            raise IOError
        return self.cifblk

    def CIFread(self,ciffile = None, cifblkname = None, cifblk = None):
        from re import sub
        from string import upper
        if ciffile != None:
            try:
                cifblk = self.CIFopen(ciffile=ciffile,cifblkname=cifblkname)
            except:
                raise IOError()
        elif cifblk == None:
            cifblk = self.cifblk

        self.atomlist.cell = [self.remove_esd(cifblk['_cell_length_a']),
                              self.remove_esd(cifblk['_cell_length_b']),
                              self.remove_esd(cifblk['_cell_length_c']),
                              self.remove_esd(cifblk['_cell_angle_alpha']),
                              self.remove_esd(cifblk['_cell_angle_beta']),
                              self.remove_esd(cifblk['_cell_angle_gamma'])]

        #self.atomlist.sgname = upper(sub("\s+","",cifblk['_symmetry_space_group_name_H-M']))
        self.atomlist.sgname = sub("\s+","",cifblk['_symmetry_space_group_name_H-M'])

        # Dispersion factors
        for i in range(len(cifblk['_atom_type_symbol'])):
            try:
                self.atomlist.dispersion[cifblk['_atom_type_symbol'][i]] =\
                    [self.remove_esd(cifblk['_atom_type_scat_dispersion_real'][i]),
                     self.remove_esd(cifblk['_atom_type_scat_dispersion_imag'][i])]
            except:
                self.atomlist.dispersion[cifblk['_atom_type_symbol'][i]] = None
                logging.warning('No dispersion factors for %s in cif file - set to zero' %cifblk['_atom_type_symbol'][i])

        for i in range(len(cifblk['_atom_site_type_symbol'])):
            label = cifblk['_atom_site_label'][i]
            #atomno = atomtype[upper(cifblk['_atom_site_type_symbol'][i])]
            atomtype = upper(cifblk['_atom_site_type_symbol'][i])
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
                logging.warning('%s has no multiplicity given in cif file >>>  Value default 1.0' %label)



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
        e.g. '1.234(56)' to '1.234'.                                                                             
        """
        from string import atof
        
        if a.find('(') == -1:
            value = atof(a)
        else:
            value = atof(a[:a.find('(')])
        return value


