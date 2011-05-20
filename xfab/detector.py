"""
The xfab.detector modules provides a numbe of function 
for calculation of detector coordinates, detector coordinate 
transformations etc.

"""

import numpy as n

def det_coor(Gt, costth, wavelength, distance, y_size, z_size, 
             dety_center, detz_center, R_tilt, tx, ty, tz):
    """
    Calculates detector coordinates dety,detz

    INPUT:
    Gt is the g-vector 
    y_size and z-size are the detector pixel size in y, and z (in microns)
    (dety_center, detz-center) is the beam center of the detector (in pixels)
    R_tilt is the rotation matrix of the detector
    (tx, ty, tz) is the position of the grain at the present omega

    OUTPUT:
    [dety, detz]
    """

    # Unit directional vector for reflection
    v = n.array([costth, 
                 wavelength/(2*n.pi)*Gt[1],
                 wavelength/(2*n.pi)*Gt[2]])
    t = (R_tilt[0, 0]*distance - \
             n.sum(R_tilt[:, 0]*n.array([tx, ty, tz])))/n.sum(R_tilt[:, 0]*v)
    Ltv = n.array([tx-distance, ty, tz])+ t*v
    dety = n.sum(R_tilt[:, 1]*Ltv)/y_size + dety_center
    detz = n.sum(R_tilt[:, 2]*Ltv)/z_size + detz_center
    return [dety, detz]

def det_coor2(tth, eta, distance, y_size, z_size, 
              dety_center, detz_center, R_tilt, tx, ty, tz,):
    """
    Calculates detector coordinates dety,detz
    Alternative to det_coor using the angles the azimuthal angle eta
    for calculation of the detector coordinates instead of the g-vector

    INPUT:
    tth is two-theta
    eta is the azimuthal angle around the beam center ccw with eta= 0 at +ylab
    y_size and z-size are the detector pixel size in y, and z (in microns)
    (dety_center, detz-center) is the beam center of the detector (in pixels)
    R_tilt is the rotation matrix of the detector
    (tx, ty, tz) is the position of the grain at the present omega

    OUTPUT:
    [dety, detz]
    """
    costth = n.cos(tth)
    sintth = n.sin(tth)
    # Unit directional vector for reflection
    v = n.array([costth, 
                 -sintth*n.sin(eta),
                 sintth*n.cos(eta)])
    t = (R_tilt[0, 0]*distance - \
             n.sum(R_tilt[:, 0]*n.array([tx, ty, tz])))/n.sum(R_tilt[:, 0]*v)
    Ltv = n.array([tx-distance, ty, tz])+ t*v
    dety = n.sum(R_tilt[:, 1]*Ltv)/y_size + dety_center
    detz = n.sum(R_tilt[:, 2]*Ltv)/z_size + detz_center
    return [dety, detz]

def det_v(Gt, costth, wavelength, distance, y_size, z_size, 
             dety_center, detz_center, R_tilt, tx, ty, tz):
    """
    Calculates direction of outgoing X-ray beam (used by fabric)
    needs to be cleaned up

    INPUT:
    Gt is the g-vector 
    y_size and z-size are the detector pixel size in y, and z (in microns)
    (dety_center, detz-center) is the beam center of the detector (in pixels)
    R_tilt is the rotation matrix of the detector
    (tx, ty, tz) is the position of the grain at the present omega

    OUTPUT:
    [dety, detz]
    """

    # Unit directional vector for reflection
    v = n.array([costth, 
                 wavelength/(2*n.pi)*Gt[1],
                 wavelength/(2*n.pi)*Gt[2]])
    return v


def detector_to_lab(dety, detz, L, py, pz, y0, z0, R_tilt):
    """
    calculates the laboratory coordinates (xl,yl,zl) from the detector
    coordinates (dety,detz), the sample-to-detector distance L, the pixel
    sizes (py and pz), the beam centre (y0,z0) and the detector tilt R_tilt 
    """
       
    lab = n.array([[L], [0], [0]]) + n.dot(R_tilt, n.array([[0],
                                                            [py*(dety-y0)],
                                                            [pz*(detz-z0)]]))
    return [lab[0][0], lab[1][0], lab[2][0]]
    

def trans_orientation(img, o11, o12, o21, o22, flipdir = 'forward'):
    """
    Transforming image matrix according to the 
    detector orientation matrix given the 
    output image  matrix will have coordinates (dety,detz)
    as defined in 
    "3DXRD and TotalCryst Geometry - Version 1.0.2" by
    H.F. Poulsen, S. Schmidt, J. Wright, H.O. Sorensen

    It is important to note that the image is not only flipped 
    using the orientation matrix, but also transposed such the 
    the returned array img has the coordinates [dety,detz].
    Therefore even the identy orientation matrix performs a transformation.
    To avoid this please use the function image_flipping
    
    Detector_orientation: [[o11,o12],[o21,o22]]
    if pretransposed to get img -> img(dety,detz) then
           [[o11,o12],[o21,o22]]
           [[  1,  0],[  0,  1]]  => nothing
           [[ -1,  0],[  0,  1]]  => fliplr
           [[  1,  0],[  0, -1]]  => flipud
           [[ -1,  0],[  0, -1]]  => flipud fliplr
    These will not be pretransposed
           [[  0,  1],[  1,  0]]  => nothing
           [[  0, -1],[ -1,  0]]  => fliplr flipud
           [[  0, -1],[  1,  0]]  => fliplr
           [[  0,  1],[ -1,  0]]  => flipud
           
    flipdir can takes the values forward or inverse
    forward: raw image -> 3DXRD standard
    inverse: 3DXRD standard -> raw image
 
    """

    if abs(o11) == 1:
        if (abs(o22) != 1) or (o12 != 0) or (o21 != 0):
            raise ValueError, 'detector orientation makes no sense 1'
        img = n.transpose(img) # to get A[i,j] be standard A[dety,detz] 
        if o11 == -1:
            if flipdir == 'forward':
                img = n.fliplr(img)
            else: #inverse direction from (dety,detz) to imageformat
                # since these direction of operations should be obversed
                # we should now do flipud before transpose 
                # But transpose(fliplr(img)) = flipud(transpose(img)) 
                img = n.flipud(img)
        if o22 == -1:
            if flipdir == 'forward':
                img = n.flipud(img)
            else: #inverse direction from (dety,detz) to imageformat
                img = n.fliplr(img)
        return img
    if abs(o12) == 1:
        if abs(o21) != 1 or (o11 != 0) or (o22 != 0):
            raise ValueError, 'detector orientation makes no sense 2'
        #transpose not needed since the matrix is transp from scratch
        if o12 == -1:
            img = n.fliplr(img)
        if o21 == -1:
            img = n.flipud(img)
        return img
    raise ValueError, 'detector orientation makes no sense 3'


def image_flipping(img, o11, o12, o21, o22, flipdir='forward'):
    """
    Transforming image matrix according to the 
    detector orientation matrix given the 
    output image  matrix will have coordinates (dety,detz)
    as defined in 
    "3DXRD and TotalCryst Geometry - Version 1.0.2" by
    H.F. Poulsen, S. Schmidt, J. Wright, H.O. Sorensen
    
    Detector_orientation: [[o11,o12],[o21,o22]]

           [[o11,o12],[o21,o22]]
           [[  1,  0],[  0,  1]]  => nothing
           [[ -1,  0],[  0,  1]]  => flipud
           [[  1,  0],[  0, -1]]  => fliplr
           [[ -1,  0],[  0, -1]]  => flipud fliplr

           [[  0,  1],[  1,  0]]  => transpose
           [[  0, -1],[ -1,  0]]  => transpose fliplr flipud
           [[  0, -1],[  1,  0]]  => transpose flipud
           [[  0,  1],[ -1,  0]]  => transpose flipud
           
    flipdir can takes the values forward or inverse
    forward: raw image -> 3DXRD standard
    inverse: 3DXRD standard -> raw image

    """

    if abs(o11) == 1:
        if (abs(o22) != 1) or (o12 != 0) or (o21 != 0):
            raise ValueError, 'detector orientation makes no sense 1'
#        img = n.transpose(img) # to get A[i,j] be standard A[dety,detz] 
        if o11 == -1:
            img = n.flipud(img)
        if o22 == -1:
            img = n.fliplr(img)
        return img
    if abs(o12) == 1:
        if abs(o21) != 1 or (o11 != 0) or (o22 != 0):
            raise ValueError, 'detector orientation makes no sense 2'
        #transpose not needed since the matrix is transp from scratch
        img = n.transpose(img) # make transpose
        
        if o12 == -1:
            if flipdir == 'forward':
                img = n.flipud(img)
            else:
                img = n.fliplr(img)
        if o21 == -1:
            if flipdir == 'forward':
                img = n.fliplr(img)
            else:
                img = n.flipud(img)
        return img
    raise ValueError, 'detector orientation makes no sense 3'



def detyz_to_xy(coor, o11, o12, o21, o22, dety_size, detz_size):
    """
    Transforming dety, detz coordinates to meet (x,y) coordinates
    of the raw image according to the  detector orientation matrix given.
    The definition of (dety,detz) is defined in 
    "3DXRD and TotalCryst Geometry - Version 1.0.2" by
    H.F. Poulsen, S. Schmidt, J. Wright, H.O. Sorensen
    
    Detector_orientation: o= [[o11,o12],
                              [o21,o22]]
    returns (dety, detz)
    """
    

    if abs(o11) == 1:
        if (abs(o22) != 1) or (o12 != 0) or (o21 != 0):
            raise ValueError, 'detector orientation makes no sense 1'
    elif abs(o12) == 1:
        if abs(o21) != 1 or (o11 != 0) or (o22 != 0):
            raise ValueError, 'detector orientation makes no sense 2'
    else:
        raise ValueError, 'detector orientation makes no sense 3'
    # transpose (dety, detz) to (detz, dety) to match order of (x,y)
    coor = n.array([coor[1], coor[0]])
    omat = n.array([[o11, o12], 
                    [o21, o22]])
    # Well we need to use the inverse operations here
    omat = n.linalg.inv(omat)
    det_size = n.array([detz_size-1, 
                        dety_size-1]) # also transpose coord in size
    coor = n.dot(omat, coor) - n.clip(n.dot(omat, det_size), 
                                      -n.max(det_size), 0)
    return coor

def xy_to_detyz(coor, o11, o12, o21, o22, dety_size, detz_size):
    """
    Transforming (x,y) coordinates of the raw image 
    to (dety, detz) coordinates according to the  
    detector orientation matrix given.
    The definition of (dety,detz) is defined in 
    "3DXRD and TotalCryst Geometry - Version 1.0.2" by
    H.F. Poulsen, S. Schmidt, J. Wright, H.O. Sorensen
    
    Detector_orientation: o= [[o11,o12],
                              [o21,o22]]
    returns (dety,detz)

    """
    

    if abs(o11) == 1:
        if (abs(o22) != 1) or (o12 != 0) or (o21 != 0):
            raise ValueError, 'detector orientation makes no sense 1'
    elif abs(o12) == 1:
        if abs(o21) != 1 or (o11 != 0) or (o22 != 0):
            raise ValueError, 'detector orientation makes no sense 2'
    else:
        raise ValueError, 'detector orientation makes no sense 3'
    omat = n.array([[o11, o12],
                    [o21, o22]])
    det_size = n.array([detz_size-1,
                        dety_size-1])
    coor = n.dot(omat, coor)- n.clip(n.dot(omat, det_size),
                                     -n.max(det_size), 0)
    # transpose (x,y) to (y,x) to match order of (dety,detz)
    coor = n.array([coor[1], coor[0]])
    return coor
     
     
def detyz_to_eta_and_radpix(coor, dety_center, detz_center):
    """
    Transforming (dety,detz) coordinates to (eta,radpix),
    where eta is in degrees (clockwise, zero point at 12 o'clock)
    and radpix is the radial number of pixels from the beam center.
    """
    
    radcoor = coor - n.array([dety_center,detz_center])
    radpix = n.sqrt(n.sum(radcoor**2))
    if radpix < 1:
        cos_eta = 1
    else:
        cos_eta = radcoor[1]/radpix
    if radcoor[0] <= 0:
        eta = 180/n.pi*n.arccos(cos_eta)
    else:
        eta = 360-180/n.pi*n.arccos(cos_eta)
        
    return [eta,radpix]
    

def eta_and_radpix_to_detyz(eta, radpix, dety_center, detz_center):
    """
    Transforming (eta,radpix) to detector coordinates (dety.detz)
    where eta is in degrees (clockwise, zero point at 12 o'clock)
    and radpix is the radial number of pixels from the beam center.
    """
    
    etarad = eta*n.pi/180.
    radcoor = radpix*n.array([-n.sin(etarad),n.cos(etarad)])
    coor = radcoor + n.array([dety_center,detz_center])
    
    return coor
    

def distort(coor, o11, o12, o21, o22, dety_size, detz_size, spatial):
    """
    To match the coordinate system of the spline file
    """

    (x, y) = detyz_to_xy(coor,
                         o11,
                         o12,
                         o21,
                         o22,
                         dety_size,
                         detz_size)

    # Do the spatial distortion
    (xd, yd) = spatial.distort(x, y)
    
    # transform coordinates back to dety,detz
    coord = xy_to_detyz([xd, yd],
                        o11,
                        o12,
                        o21,
                        o22,
                        dety_size,
                        detz_size)
    return coord
