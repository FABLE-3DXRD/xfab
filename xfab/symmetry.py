"""
xfab.symmetry has a few function for calculation of symmetry
related stuff
"""

from numpy import zeros, arccos, pi, dot, transpose, array, concatenate, cos, sin, dot
from numpy.linalg import det, inv
from xfab import tools


def Umis(umat_1, umat_2, crystal_system):
    """    
     Determines misorientations between grain orientations     
     Input:    umat_1 and umat_2 orientation matrices   
               crystal_system  1: Triclinic
                               2: Monoclinic
                               3: Orthorhombic
                               4: Tetragonal
                               5: Trigonal
                               6: Hexagonal
                               7: Cubic
     Output:   (rotation number, misorientation in deg)
    """

    #rotation matrix defined by crystal system 
    rot = rotations(crystal_system)
    nrot = rot.shape[0]
    t_save = zeros((0, 2))
    for k in range(nrot):
        g_vector = dot(umat_2, transpose(dot(umat_1, rot[k])))
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

    
def add_perm(hkl, crystal_system):
    """
    apply the permutation according to the crystal system on the 
    hkl and print the result
    """
    
    perm = permutations(crystal_system)
    nperm = perm.shape[0]
    
    for k in range(nperm):
        print dot(perm[k],hkl)

        
def add_rot(umat, crystal_system):
    """
    apply the rotation according to the crystal system on the 
    orientation matrix U and print the result
    """
    
    rot = rotations(crystal_system)
    nrot = rot.shape[0]
    
    for k in range(nrot):
        print dot(umat, rot[k])

        
def permutations(crystal_system):
    
    """ 
    permutations returns the set of indistinguasible lattice permutations
    
    perm = permutations(crystal_system)
    hkl_new = perm*hkl
    
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

    
def rotations(crystal_system):
    
    """ 
    rotations returns the set of unitary rotation matrices
    corresponding to the indistinguasible lattice permutations
    The values of the function only differ from permutations for 
    trigonal and hexagonal crystal systems
    
    rot = rotations(crystal_system)
    U_new = U*rot
    
    U_new*B*perm*hkl = U*B*hkl, so U_new = U*B*perminv*Binv and
    rot = B*perminv*Binv. If B or perm is diagonal, rot = perminv.
    perminv included in perm, but to get the i'th entry of rot
    to correspond to the i'th entry of perm one must use perminv.
    
    crystal_system can be one of the following values   
    1: Triclinic
    2: Monoclinic
    3: Orthorhombic
    4: Tetragonal
    5: Trigonal
    6: Hexagonal
    7: Cubic
    
    Jette Oddershede, Riso, 8/2/2010
    """

    if crystal_system < 1 or crystal_system > 7:
        raise ValueError('Crystal system shoud have a value between 1 and 7')

    if crystal_system == 1: # Triclinic
        rot = permutations(crystal_system)
        for i in range(len(rot)):       
            rot[i] = inv(rot[i])

    if crystal_system == 2: # Monoclinic
        rot = permutations(crystal_system)
        for i in range(len(rot)):       
            rot[i] = inv(rot[i])

    if crystal_system == 3: # Orthorhombic
        rot = permutations(crystal_system)
        for i in range(len(rot)):       
            rot[i] = inv(rot[i])
 
    if crystal_system == 4: # Tetragonal
        rot = permutations(crystal_system)
        for i in range(len(rot)):       
            rot[i] = inv(rot[i])

    if crystal_system == 5: # Trigonal
        perm = permutations(crystal_system)
        B = tools.form_b_mat([1.,1.,1.,90.,90.,120.])
        Binv = inv(B)
        rot = zeros((6, 3, 3))
        for i in range(len(perm)):       
            rot[i] = dot(B,dot(inv(perm[i]),Binv))
                
    if crystal_system == 6: # Hexagonal
        perm = permutations(crystal_system)
        B = tools.form_b_mat([1.,1.,1.,90.,90.,120.])
        Binv = inv(B)
        rot = zeros((12, 3, 3))
        for i in range(len(perm)):       
            rot[i] = dot(B,dot(inv(perm[i]),Binv))

    if crystal_system == 7: # Cubic
        rot = permutations(crystal_system)
        for i in range(len(rot)):       
            rot[i] = inv(rot[i])

    return rot
