
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
