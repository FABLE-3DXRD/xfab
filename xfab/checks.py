
"""The checks module implements input checks for various xfab functions. These checks
are meant to be called upon during runtime to catch errors made on the user side. To
toogle checks in xfab a package wide constant CHECKS defined in the __init__.py should
be accessed as:

import xfab
xfab.CHECKS.activated = False

"""
import numpy as np

def _check_rotation_matrix(U):
    """Verify that a 3 x 3 matrix is a rotation matrix.

    Args:
        U: 3x3 matrix represented as a numpy array of shape=(3,3).

    """
    if not np.allclose( np.dot(U.T, U), np.eye(3,3), atol=1e-6):
        raise ValueError("orientation matrix U is not unitary, np.dot(U.T, U)!=np.eye(3,3)")

    if not np.allclose( np.linalg.det(U), 1.0 ):
        raise ValueError("orientation matrix U has a non unity determinant np.linalg.det(U)!=1.0")

def _check_euler_angles(phi1, PHI, phi2):
    """Verify that all three Euler angles lies in the range [0,2*pi].

    Args:
        phi1, PHI, and phi2: Euler angles in radians.

    """
    if not (0<=phi1<=np.pi*2):
        raise ValueError("Euler angle phi1="+str(phi1)+" is not in range [0,2*pi]")

    if not (0<=PHI<=np.pi*2):
        raise ValueError("Euler angle PHI="+str(PHI)+" is not in range [0,2*pi]")

    if not (0<=phi2<=np.pi*2):
        raise ValueError("Euler angle phi2="+str(phi2)+" is not in range [0,2*pi]")

class _checkState(object):

    """A checkState object can be used to determine if checks are active in xfab
        The default mode is active and can be togled. This allows users to in some
        specific scenarios speed up their codes. Furthermore, if __debug__==False
        all checks are inactivated, i.e using python -O to run codes will remove 
        any checks.
    """
    def __init__(self):
        self._run_checks = True

    @property
    def activated(self):
        """True if checks are to be run.
        """
        return self._run_checks and __debug__

    @activated.setter
    def activated(self, value):
        if value is not True and value is not False:
            raise ValueError("Please supply a boolean True or False")
        else:
            self._run_checks = value