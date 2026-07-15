"""
Validate form factor database against expected values.
f(0) = sum(a_i) + c should equal Z (atomic number) for neutral atoms.
"""

import pytest
from xfab import atomlib

Z = {
    'H': 1,    'HE': 2,   'LI': 3,   'BE': 4,   'B': 5,
    'C': 6,    'N': 7,    'O': 8,    'F': 9,    'NE': 10,
    'NA': 11,  'MG': 12,  'AL': 13,  'SI': 14,  'P': 15,
    'S': 16,   'CL': 17,  'AR': 18,  'K': 19,   'CA': 20,
    'SC': 21,  'TI': 22,  'V': 23,   'CR': 24,  'MN': 25,
    'FE': 26,  'CO': 27,  'NI': 28,  'CU': 29,  'ZN': 30,
    'GA': 31,  'GE': 32,  'AS': 33,  'SE': 34,  'BR': 35,
    'KR': 36,  'RB': 37,  'SR': 38,  'Y': 39,   'ZR': 40,
    'NB': 41,  'MO': 42,  'TC': 43,  'RU': 44,  'RH': 45,
    'PD': 46,  'AG': 47,  'CD': 48,  'IN': 49,  'SN': 50,
    'SB': 51,  'TE': 52,  'I': 53,   'XE': 54,  'CS': 55,
    'BA': 56,  'LA': 57,  'CE': 58,  'PR': 59,  'ND': 60,
    'PM': 61,  'SM': 62,  'EU': 63,  'GD': 64,  'TB': 65,
    'DY': 66,  'HO': 67,  'ER': 68,  'TM': 69,  'YB': 70,
    'LU': 71,  'HF': 72,  'TA': 73,  'W': 74,   'RE': 75,
    'OS': 76,  'IR': 77,  'PT': 78,  'AU': 79,  'HG': 80,
    'TL': 81,  'PB': 82,  'BI': 83,  'PO': 84,  'AT': 85,
    'RN': 86,  'FR': 87,  'RA': 88,  'AC': 89,  'TH': 90,
    'PA': 91,  'U': 92,   'NP': 93,  'PU': 94,
}


@pytest.mark.parametrize("key", sorted(Z.keys()))
def test_form_factor_f0_equals_Z(key):
    coeffs = atomlib.formfactor[key]
    a_sum = coeffs[0] + coeffs[1] + coeffs[2] + coeffs[3]
    c = coeffs[8]
    f0 = a_sum + c
    assert abs(f0 - Z[key]) < 0.5, \
        "{}: f(0) = {:.3f}, expected Z = {}, diff = {:.3f}".format(
            key, f0, Z[key], f0 - Z[key])


def test_all_elements_present():
    assert len(atomlib.formfactor) == len(Z), \
        "Expected {} elements, got {}".format(len(Z), len(atomlib.formfactor))
    for key in Z:
        assert key in atomlib.formfactor, "Missing element: {}".format(key)


def test_form_factor_decreases_with_stl():
    import numpy as np
    from xfab.structure import FormFactor

    stl_vals = [0.0, 0.2, 0.5, 1.0]
    for key in ['PB', 'TI', 'O', 'FE', 'AU']:
        prev = None
        for stl in stl_vals:
            f = FormFactor(key, stl)
            if prev is not None:
                assert f < prev, \
                    "{}: f({})={:.3f} >= f(prev_stl)={:.3f}".format(
                        key, stl, f, prev)
            prev = f
