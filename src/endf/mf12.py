# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO
from warnings import warn

from .records import get_head_record, get_tab1_record, get_list_record


def parse_mf12(file_obj: TextIO) -> dict:
    """Parse photon production multiplicities from MF=12

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Photon production multiplicity / transition probability data

    """
    ZA, AWR, LO, LG, NK, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LO': LO, 'NK': NK}

    if LO == 1:
        # Multiplicities given -- start by reading total yield
        if NK > 1:
            _, data['Y'] = get_tab1_record(file_obj)

        # Read multiplicities
        data['multiplicities'] = []
        for k in range(NK):
            (Eg, ES, LP, LF), y = get_tab1_record(file_obj)
            data_k = {'Eg': Eg, 'ES': ES, 'LP': LP, 'LF': LF, 'y': y}
            data['multiplicities'].append(data_k)

    elif LO == 2:
        # Store whether simple (LG=1) or complex (LG=2) transitions
        data['LG'] = LG

        # Get transition probability data
        (ES_NS, _, LP, _, _, NT), values = get_list_record(file_obj)
        data['ES_NS'] = ES_NS
        data['LP'] = LP
        data['NT'] = NT
        data['transitions'] = transition = []
        for i in range(NT):
            if LG == 1:
                ES, TP = values[2*i : 2*(i + 1)]
                transition.append({'ES': ES, 'TP': TP})
            elif LG == 2:
                ES, TP, GP = values[3*i : 3*(i + 1)]
                transition.append({'ES': ES, 'TP': TP, 'GP': GP})
    else:
        warn(f"Unrecognized LO value: {LO}")

    return data
