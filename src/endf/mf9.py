# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf9_mf10(file_obj: TextIO, MF: int) -> dict:
    """Parse radionuclide production data from MF=9 or MF=10

    Parameters
    ----------
    file_obj
        File-like object to read from
    MF
        File number

    Returns
    -------
    dict
        Radionuclide production data

    """
    ZA, AWR, LIS, _, NS, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LIS': LIS, 'NS': NS}
    data['levels'] = []
    for _ in range(NS):
        # Determine what the product is
        (QM, QI, IZAP, LFS), func = get_tab1_record(file_obj)
        level_data = {'QM': QM, 'QI': QI, 'IZAP': IZAP, 'LFS': LFS}
        if MF == 9:
            level_data['Y'] = func
        else:
            level_data['sigma'] = func
        data['levels'].append(level_data)

    return data
