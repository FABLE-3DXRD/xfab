"""
xfab.symmetry has a few function for calculation of symmetry
related stuff
"""

from numpy import zeros, arccos, pi, dot, transpose, array, concatenate
from numpy.linalg import det


def Umis(umat_1, umat_2, crystal_system):
    """
    
     Determines misorientations between grain orientations 
    
     calc(ubifile,crystal_system)
    
     Input:    graininfo.mat is the simul_farfield file containing simulated
               orientations
    
               ubifile: file with orientations as (UB)^{-1} from
                        GrainSpotter or ImageD11-index
    
               crystal_system  1: Triclinic
                               2: Monoclinic
                               3: Orthorhombic
                               4: Tetragonal
                               5: Trigonal
                               6: Hexagonal
                               7: Cubic

    """

    #Permutation matrix defined by crystal system 
    perm = permutations(crystal_system)
    nperm = perm.shape[0]
    t_save = zeros((0, 2))
    for k in range(nperm):
        g_vector = dot(umat_2, transpose(dot(umat_1, perm[k])))
        detg = det(g_vector)
        if detg < 0.9999 or detg > 1.0001:
            print 'mistake %f' % detg
        else:
            length = (g_vector.diagonal().sum()-1.)/2.
            if abs(length) > 1.00000000:
                if length > 1:
                    length = 1.
                else:
                    length = -1.                    
            misangle = arccos(length)
            t_save = concatenate((t_save, array([[k, misangle*180/pi]])), 0)
    return t_save

    
def add_perm(umat, crystal_system):
    """
    apply the permutation according to the crystal system on the 
    orientation matrix U and print the result
    """
    perm = permutations(crystal_system)
    nperm = perm.shape[0]
    
    for k in range(nperm):
        print dot(umat, perm[k])

def permutations(crystal_system):
    
    """ 
    lattperm returns the set of indistinguasible lattice permutations
    
    [perm] = lattperm(crystal_system)
    
    crystal_system can be one of the following values
    
    1: Triclinic
    2: Monoclinic
    3: Orthorhombic
    4: Tetragonal
    5: Trigonal
    6: Hexagonal
    7: Cubic
    
    Henning Osholm Sorensen, Riso, 30/6/2006
    Implemented in python, 12/7/2008
    """

    if crystal_system < 1 or crystal_system > 7:
        raise ValueError('Crystal system shoud have a value between 1 and 7')

    if crystal_system == 1: # Triclinic
        perm = zeros((1, 3, 3))
        perm[0]  = [[ 1., 0, 0], [ 0, 1., 0], [ 0, 0, 1.]]

    if crystal_system == 2: # Monoclinic
        perm = zeros((2, 3, 3))
        perm[0]  = [[ 1, 0, 0], [ 0, 1, 0], [ 0, 0,  1]]
        perm[1]  = [[-1, 0, 0], [ 0, 1, 0], [ 0, 0, -1]]

    if crystal_system == 3: # Orthorhombic
        perm = zeros((4, 3, 3))
        perm[0]  = [[ 1, 0, 0], [ 0,  1, 0], [ 0, 0,  1]]
        perm[1]  = [[-1, 0, 0], [ 0, -1, 0], [ 0, 0,  1]]
        perm[2]  = [[-1, 0, 0], [ 0,  1, 0], [ 0, 0, -1]]
        perm[3]  = [[ 1, 0, 0], [ 0, -1, 0], [ 0, 0, -1]]
 
    if crystal_system == 4: # Tetragonal
        perm = zeros((8, 3, 3))
        perm[0]  = [[ 1,  0, 0], [ 0,  1, 0], [ 0, 0,  1]]
        perm[1]  = [[-1,  0, 0], [ 0, -1, 0], [ 0, 0,  1]]
        perm[2]  = [[ 0,  1, 0], [-1,  0, 0], [ 0, 0,  1]]
        perm[3]  = [[ 0, -1, 0], [ 1,  0, 0], [ 0, 0,  1]]
        perm[4]  = [[-1,  0, 0], [ 0,  1, 0], [ 0, 0, -1]]
        perm[5]  = [[ 1,  0, 0], [ 0, -1, 0], [ 0, 0, -1]]
        perm[6]  = [[ 0,  1, 0], [ 1,  0, 0], [ 0, 0, -1]]
        perm[7]  = [[ 0, -1, 0], [-1,  0, 0], [ 0, 0, -1]]

    if crystal_system == 5: # Trigonal
        perm = zeros((6, 3, 3))
        perm[0]  = [[ 1,  0, 0], [ 0,  1, 0], [ 0, 0,  1]]
        perm[1]  = [[ 0,  1, 0], [-1, -1, 0], [ 0, 0,  1]]
        perm[2]  = [[-1, -1, 0], [ 1,  0, 0], [ 0, 0,  1]]
        perm[3]  = [[ 0,  1, 0], [ 1,  0, 0], [ 0, 0, -1]]
        perm[4]  = [[ 1,  0, 0], [-1, -1, 0], [ 0, 0, -1]]
        perm[5]  = [[-1, -1, 0], [ 0,  1, 0], [ 0, 0, -1]]

    if crystal_system == 6: # Hexagonal
        perm = zeros((12, 3, 3))
        perm[0]  = [[ 1,  0, 0], [ 0,  1, 0], [ 0, 0,  1]]
        perm[1]  = [[ 0,  1, 0], [-1, -1, 0], [ 0, 0,  1]]
        perm[2]  = [[-1, -1, 0], [ 1,  0, 0], [ 0, 0,  1]]
        perm[3]  = [[-1,  0, 0], [ 0, -1, 0], [ 0, 0,  1]]
        perm[4]  = [[ 0, -1, 0], [ 1,  1, 0], [ 0, 0,  1]]
        perm[5]  = [[ 1,  1, 0], [-1,  0, 0], [ 0, 0,  1]]
        perm[6]  = [[ 0,  1, 0], [ 1,  0, 0], [ 0, 0, -1]]
        perm[7]  = [[ 1,  0, 0], [-1, -1, 0], [ 0, 0, -1]]
        perm[8]  = [[-1, -1, 0], [ 0,  1, 0], [ 0, 0, -1]]
        perm[9]  = [[ 0, -1, 0], [-1,  0, 0], [ 0, 0, -1]]
        perm[10] = [[-1,  0, 0], [ 1,  1, 0], [ 0, 0, -1]]
        perm[11] = [[ 1,  1, 0], [ 0, -1, 0], [ 0, 0, -1]]

    if crystal_system == 7: # Cubic
        perm = zeros((24, 3, 3))
        perm[0]  = [[ 1,  0,  0], [ 0,  1,  0], [ 0,  0,  1]]
        perm[1]  = [[ 1,  0,  0], [ 0, -1,  0], [ 0,  0, -1]]
        perm[2]  = [[ 1,  0,  0], [ 0,  0, -1], [ 0,  1,  0]]
        perm[3]  = [[ 1,  0,  0], [ 0,  0,  1], [ 0, -1,  0]]
        perm[4]  = [[-1,  0,  0], [ 0,  1,  0], [ 0,  0, -1]]
        perm[5]  = [[-1,  0,  0], [ 0, -1,  0], [ 0,  0,  1]]
        perm[6]  = [[-1,  0,  0], [ 0,  0, -1], [ 0, -1,  0]]
        perm[7]  = [[-1,  0,  0], [ 0,  0,  1], [ 0,  1,  0]]
        perm[8]  = [[ 0,  1,  0], [-1,  0,  0], [ 0,  0,  1]]
        perm[9]  = [[ 0,  1,  0], [ 0,  0, -1], [-1,  0,  0]]
        perm[10] = [[ 0,  1,  0], [ 1,  0,  0], [ 0,  0, -1]]
        perm[11] = [[ 0,  1,  0], [ 0,  0,  1], [ 1,  0,  0]]
        perm[12] = [[ 0, -1,  0], [ 1,  0,  0], [ 0,  0,  1]]
        perm[13] = [[ 0, -1,  0], [ 0,  0, -1], [ 1,  0,  0]]
        perm[14] = [[ 0, -1,  0], [-1,  0,  0], [ 0,  0, -1]]
        perm[15] = [[ 0, -1,  0], [ 0,  0,  1], [-1,  0,  0]]
        perm[16] = [[ 0,  0,  1], [ 0,  1,  0], [-1,  0,  0]]
        perm[17] = [[ 0,  0,  1], [ 1,  0,  0], [ 0,  1,  0]]
        perm[18] = [[ 0,  0,  1], [ 0, -1,  0], [ 1,  0,  0]]
        perm[19] = [[ 0,  0,  1], [-1,  0,  0], [ 0, -1,  0]]
        perm[20] = [[ 0,  0, -1], [ 0,  1,  0], [ 1,  0,  0]]
        perm[21] = [[ 0,  0, -1], [-1,  0,  0], [ 0,  1,  0]]
        perm[22] = [[ 0,  0, -1], [ 0, -1,  0], [-1,  0,  0]]
        perm[23] = [[ 0,  0, -1], [ 1,  0,  0], [ 0, -1,  0]]

    return perm
