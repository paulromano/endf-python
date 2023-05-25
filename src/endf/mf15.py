# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_head_record, get_tab1_record, get_tab2_record


def parse_mf15(file_obj: TextIO) -> dict:
    """Parse continuous photon energy spectra from MF=15

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Continuous photon energy spectra data

    """
    ZA, AWR, _, _, NC, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'NC': NC}
    data['subsections'] = []
    for j in range(NC):
        # Read probability for j-th partial distribution
        params, p = get_tab1_record(file_obj)
        LF = params[3]
        subsec = {'LF': LF, 'p': p}

        # Read tabulated distributions
        (*_, NE), subsec['E_int'] = get_tab2_record(file_obj)
        subsec['NE'] = NE
        subsec['E'] = np.empty(NE)
        subsec['g'] = []
        for i in range(NE):
            (_, E, _, _), g = get_tab1_record(file_obj)
            subsec['E'][i] = E
            subsec['g'].append(g)
        data['subsections'].append(subsec)

    return data
