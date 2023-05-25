# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_head_record, get_cont_record, get_tab2_record, \
    get_list_record, get_tab1_record


def parse_mf14(file_obj: TextIO) -> dict:
    """Parse photon angular distributions from MF=14

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Photon angular distribution data

    """
    ZA, AWR, LI, LTT, NK, NI = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LI': LI, 'NK': NK}

    # If all photons are isotropic, exit early
    if LI == 1:
        return data

    data['LTT'] = LTT
    data['NI'] = NI
    data['subsections'] = []

    # Subsections for isotropic photons
    for i in range(NI):
        EG, ES, *_ = get_cont_record(file_obj)
        data['subsections'].append({'EG': EG, 'ES': ES})

    # Subsections for anisotropic photons
    for i in range(NI, NK):
        (EG, ES, _, _, NR, NE), E_int = get_tab2_record(file_obj)
        subsec = {'EG': EG, 'ES': ES, 'NE': NE, 'E_int': E_int}
        subsec['E'] = np.empty(NE)
        if LTT == 1:
            # Legendre coefficient representation
            subsec['NL'] = np.empty(NE)
            subsec['a_lk'] = []
            for i in range(NE):
                (_, E, _, _, NL, _ ), a_lk = get_list_record(file_obj)
                subsec['E'][i] = E
                subsec['NL'][i] = NL
                subsec['a_lk'].append(a_lk)
        elif LTT == 2:
            # Tabulated angular distribution
            subsec['p_k'] = []
            for i in range(NE):
                (_, E, *_), p_k = get_tab1_record(file_obj)
                subsec['E'][i] = E
                subsec['p_k'].append(p_k)

        data['subsections'].append(subsec)

    return data
