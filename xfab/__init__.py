"""
xfab is a library of data and functions for use in Crystallographic computations. xfab is part of the project fable: http://fable.wiki.sourceforge.net
"""


"""Set a package wide global variable that toogle input/output on functions.
Checks are on by default. Toogling of checks can be achived as:

    import xfab
    xfab.CHECKS.activated = False

"""
from xfab.checks import _checkState
CHECKS = _checkState()
