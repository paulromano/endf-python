# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf27(file_obj: TextIO) -> dict:
    """Parse atomic form factors / scattering functions from MF=27

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Atomic form factor or scattering function data

    """
    ZA, AWR, *_ = get_head_record(file_obj)
    params, H = get_tab1_record(file_obj)
    return {'ZA': ZA, 'AWR': AWR, 'Z': params[1], 'H': H}
