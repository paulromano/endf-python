# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf23(file_obj: TextIO) -> dict:
    """Parse photon cross sections from MF=23

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Photon cross section data

    """
    ZA, AWR, *_ = get_head_record(file_obj)
    params, xs = get_tab1_record(file_obj)
    return {
        'ZA': ZA,
        'AWR': AWR,
        'EPE': params[0],
        'EFL': params[1],
        'sigma': xs
    }
