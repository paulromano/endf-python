# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .mf33 import parse_mf33_subsection
from .records import get_head_record, get_cont_record


def parse_mf40(file_obj: TextIO) -> dict:
    """Parse covariances of radionuclide production from MF=40

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Radionuclide production covariance data

    """
    ZA, AWR, LIS, _, NS, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LIS': LIS, 'NS': NS, 'subsections': []}
    for _ in range(NS):
        QM, QI, IZAP, LFS, _, NL = get_cont_record(file_obj)
        subsection = {'QM': QM, 'QI': QI, 'IZAP': IZAP, 'LFS': LFS, 'NL': NL}
        subsection['subsubsections'] = []
        for _ in range(NL):
            # Each sub-subsection has same format as in MF=33
            subsubsection = parse_mf33_subsection(file_obj)
            subsection['subsubsections'].append(subsubsection)

        data['subsections'].append(subsection)

    return data
