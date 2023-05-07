import re
from typing import TextIO, Tuple

import numpy as np

from .function import Tabulated1D, Tabulated2D
from ._records import float_endf

ENDF_FLOAT_RE = re.compile(r'([\s\-\+]?\d*\.\d+)([\+\-]) ?(\d+)')


def py_float_endf(s: str) -> float:
    """Convert string of floating point number in ENDF to float.

    The ENDF-6 format uses an 'e-less' floating point number format,
    e.g. -1.23481+10. Trying to convert using the float built-in won't work
    because of the lack of an 'e'. This function allows such strings to be
    converted while still allowing numbers that are not in exponential notation
    to be converted as well.

    Parameters
    ----------
    s : str
        Floating-point number from an ENDF file

    Returns
    -------
    float
        The number

    """
    return float(ENDF_FLOAT_RE.sub(r'\1e\2\3', s))


def int_endf(s: str) -> int:
    """Convert string of integer number in ENDF to int.

    The ENDF-6 format technically allows integers to be represented by a field
    of all blanks. This function acts like int(s) except when s is a string of
    all whitespace, in which case zero is returned.

    Parameters
    ----------
    s : str
        Integer or spaces

    Returns
    -------
    integer
        The number or 0
    """
    return 0 if s.isspace() else int(s)


def get_text_record(file_obj) -> str:
    """Return data from a TEXT record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    str
        Text within the TEXT record

    """
    return file_obj.readline()[:66]


def get_cont_record(file_obj, skip_c=False):
    """Return data from a CONT record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from
    skip_c : bool
        Determine whether to skip the first two quantities (C1, C2) of the CONT
        record.

    Returns
    -------
    tuple
        The six items within the CONT record

    """
    line = file_obj.readline()
    if skip_c:
        C1 = None
        C2 = None
    else:
        C1 = float_endf(line[:11])
        C2 = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    N1 = int_endf(line[44:55])
    N2 = int_endf(line[55:66])
    return (C1, C2, L1, L2, N1, N2)


def get_head_record(file_obj):
    """Return data from a HEAD record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    tuple
        The six items within the HEAD record

    """
    line = file_obj.readline()
    ZA = int(float_endf(line[:11]))
    AWR = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    N1 = int_endf(line[44:55])
    N2 = int_endf(line[55:66])
    return (ZA, AWR, L1, L2, N1, N2)


def get_list_record(file_obj: TextIO) -> Tuple[list, np.ndarray]:
    """Return data from a LIST record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    list
        The six items within the header
    numpy.ndarray
        The values within the list

    """
    # determine how many items are in list
    items = get_cont_record(file_obj)
    NPL = items[4]

    # read items
    b = np.empty(NPL)
    for i in range((NPL - 1)//6 + 1):
        line = file_obj.readline()
        n = min(6, NPL - 6*i)
        for j in range(n):
            b[6*i + j] = float_endf(line[11*j:11*(j + 1)])

    return (items, b)


def get_tab1_record(file_obj):
    """Return data from a TAB1 record in an ENDF-6 file.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    list
        The six items within the header
    openmc.data.Tabulated1D
        The tabulated function

    """
    # Determine how many interpolation regions and total points there are
    line = file_obj.readline()
    C1 = float_endf(line[:11])
    C2 = float_endf(line[11:22])
    L1 = int_endf(line[22:33])
    L2 = int_endf(line[33:44])
    n_regions = int_endf(line[44:55])
    n_pairs = int_endf(line[55:66])
    params = [C1, C2, L1, L2]

    # Read the interpolation region data, namely NBT and INT
    breakpoints = np.zeros(n_regions, dtype=int)
    interpolation = np.zeros(n_regions, dtype=int)
    m = 0
    for i in range((n_regions - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_regions - m)
        for j in range(to_read):
            breakpoints[m] = int_endf(line[0:11])
            interpolation[m] = int_endf(line[11:22])
            line = line[22:]
            m += 1

    # Read tabulated pairs x(n) and y(n)
    x = np.zeros(n_pairs)
    y = np.zeros(n_pairs)
    m = 0
    for i in range((n_pairs - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_pairs - m)
        for j in range(to_read):
            x[m] = float_endf(line[:11])
            y[m] = float_endf(line[11:22])
            line = line[22:]
            m += 1

    return params, Tabulated1D(x, y, breakpoints, interpolation)


def get_tab2_record(file_obj):
    # Determine how many interpolation regions and total points there are
    params = get_cont_record(file_obj)
    n_regions = params[4]

    # Read the interpolation region data, namely NBT and INT
    breakpoints = np.zeros(n_regions, dtype=int)
    interpolation = np.zeros(n_regions, dtype=int)
    m = 0
    for _ in range((n_regions - 1)//3 + 1):
        line = file_obj.readline()
        to_read = min(3, n_regions - m)
        for _ in range(to_read):
            breakpoints[m] = int(line[0:11])
            interpolation[m] = int(line[11:22])
            line = line[22:]
            m += 1

    return params, Tabulated2D(breakpoints, interpolation)


def get_intg_record(file_obj):
    """
    Return data from an INTG record in an ENDF-6 file. Used to store the
    covariance matrix in a compact format.

    Parameters
    ----------
    file_obj : file-like object
        ENDF-6 file to read from

    Returns
    -------
    numpy.ndarray
        The correlation matrix described in the INTG record
    """
    # determine how many items are in list and NDIGIT
    items = get_cont_record(file_obj)
    ndigit = items[2]
    npar = items[3]    # Number of parameters
    nlines = items[4]  # Lines to read
    NROW_RULES = {2: 18, 3: 12, 4: 11, 5: 9, 6: 8}
    nrow = NROW_RULES[ndigit]

    # read lines and build correlation matrix
    corr = np.identity(npar)
    for i in range(nlines):
        line = file_obj.readline()
        ii = int_endf(line[:5]) - 1  # -1 to account for 0 indexing
        jj = int_endf(line[5:10]) - 1
        factor = 10**ndigit
        for j in range(nrow):
            if jj+j >= ii:
                break
            element = int_endf(line[11+(ndigit+1)*j:11+(ndigit+1)*(j+1)])
            if element > 0:
                corr[ii, jj] = (element+0.5)/factor
            elif element < 0:
                corr[ii, jj] = (element-0.5)/factor

    # Symmetrize the correlation matrix
    corr = corr + corr.T - np.diag(corr.diagonal())
    return corr
