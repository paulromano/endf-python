# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf3(file_obj: TextIO) -> dict:
    """Parse reaction cross sections from MF=3

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Cross section data

    """
    ZA, AWR, *_ = get_head_record(file_obj)
    params, xs = get_tab1_record(file_obj)
    return {
        'ZA': ZA,
        'AWR': AWR,
        'QM': params[0],
        'QI': params[1],
        'LR': params[3],
        'sigma': xs
    }
