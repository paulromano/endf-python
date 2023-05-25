# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_head_record, get_cont_record, get_list_record


def parse_mf34(file_obj: TextIO, MT: int) -> dict:
    """Parse covariances of angular distributions from MF=34

    Parameters
    ----------
    file_obj
        File-like object to read from
    MT
        Reaction number

    Returns
    -------
    dict
        Angular distribution covariance data

    """
    ZA, AWR, _, LTT, _, NMT1 = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LTT': LTT, 'NMT1': NMT1, 'subsections': []}
    for _ in range(NMT1):
        _, _, MAT1, MT1, NL, NL1 = get_cont_record(file_obj)
        if MT1 == 0 or MT == MT1:
            NSS = NL*(NL + 1)//2
        else:
            NSS = NL*NL1
        subsection = {'MAT1': MAT1, 'MT1': MT1, 'NL': NL, 'NSS': NSS}
        subsection['L'] = np.empty(NSS)
        subsection['L1'] = np.empty(NSS)
        subsection['NI'] = np.empty(NSS)
        subsection['subsubsections'] = []

        for n in range(NSS):
            _, _, L, L1, LCT, NI = get_cont_record(file_obj)
            subsection['L'][n] = L
            subsection['L1'][n] = L1
            subsection['NI'][n] = NI
            if n == 0:
                subsection['LCT'] = LCT

            subsub = {
                'LS': np.empty(NI),
                'LB': np.empty(NI),
                'NT': np.empty(NI),
                'NE': np.empty(NI),
                'Data': []
            }
            for m in range(NI):
                (_, _, LS, LB, NT, NE), values = get_list_record(file_obj)
                subsub['LS'][m] = LS
                subsub['LB'][m] = LS
                subsub['NT'][m] = NT
                subsub['NE'][m] = NE
                subsub['Data'].append(values)
            subsection['subsubsections'].append(subsub)

    return data
