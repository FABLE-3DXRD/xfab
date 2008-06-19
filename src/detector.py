
import numpy as n


def det_coor(Gt, costth, wavelength, distance, y_size, z_size, dety_center, detz_center,
             R_tilt, tx, ty, tz,):
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

    v = n.array([costth, 
                 wavelength/(2*n.pi)*Gt[1],
                 wavelength/(2*n.pi)*Gt[2]])
    t = R_tilt[0,0]*distance/n.sum(R_tilt[:,0]*v)
    Ltv = n.array([tx-distance, ty, tz])+ t*v
    dety = n.sum(R_tilt[:,1]*Ltv)/y_size + dety_center
    detz = n.sum(R_tilt[:,2]*Ltv)/z_size + detz_center
    return [dety, detz]


def trans_orientation(img,o11,o12,o21,o22,dir='forward'):
    """
    Transforming image matrix according to the 
    detector orientation matrix given the 
    output image  matrix will have coordinates (dety,detz)
    as defined in 
    "3DXRD and TotalCryst Geometry - Version 1.0.2" by
    H.F. Poulsen, S. Schmidt, J. Wright, H.O. Sorensen
    
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
    """

    if abs(o11) == 1:
        if (abs(o22) != 1) or (o12 != 0) or (o21 != 0):
            raise ValueError, 'detector orientation makes no sense 1'
        img = n.transpose(img) # to get A[i,j] be standard A[dety,detz] 
        if o11 == -1:
            if dir == 'forward':
                img = n.fliplr(img)
            else: #inverse direction from (dety,detz) to imageformat
                # since these direction of operations should be obversed
                # we should now do flipud before transpose 
                # But transpose(fliplr(img)) = flipud(transpose(img)) 
                img = n.flipud(img)
        if o22 == -1:
            if dir == 'forward':
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



def detyz2xy(coor,o11,o12,o21,o22,dety_size,detz_size):
    """
    Transforming dety, detz coordinates to meet (x,y) coordinates
    of the raw image according to the  detector orientation matrix given.
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
    # transpose (dety,detz) to (detz,dety) to match order of (x,y)
    coor = n.array([coor[1],coor[0]])
    o = n.array([[o11,o12],[o21,o22]])
    det_size = n.array([dety_size-1,detz_size-1])
    coor = n.dot(o,coor)- n.clip(n.dot(o,det_size),-n.max(det_size),0)
    return coor

def xy2detyz(coor,o11,o12,o21,o22,dety_size,detz_size):
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
    o = n.array([[o11,o12],[o21,o22]])
    det_size = n.array([dety_size-1,detz_size-1])
    coor = n.dot(o,coor)- n.clip(n.dot(o,det_size),-n.max(det_size),0)
    # transpose (x,y) to (y,x) to match order of (dety,detz)
    coor = n.array([coor[1],coor[0]])
    return coor
