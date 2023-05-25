# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf13(file_obj: TextIO) -> dict:
    """Parse photon production cross sections from MF=13

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Photon production cross section data

    """
    ZA, AWR, _, _, NK, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'NK': NK}

    # Read total photon production
    if NK > 1:
        _, data['sigma_total'] = get_tab1_record(file_obj)

    # Read production cross sections for each photon
    data['photons'] = []
    for k in range(NK):
        (EG, ES, LP, LF), sigma = get_tab1_record(file_obj)
        photon = {'EG': EG, 'ES': ES, 'LP': LP, 'LF': LF, 'sigma': sigma}
        data['photons'].append(photon)

    return data
