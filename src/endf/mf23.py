# SPDX-FileCopyrightText: 2023 International Atomic Energy Agency
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf23(file_obj: TextIO) -> dict:
    # Generate cross section
    ZA, AWR, *_ = get_head_record(file_obj)
    params, xs = get_tab1_record(file_obj)
    return {
        'ZA': ZA,
        'AWR': AWR,
        'EPE': params[0],
        'EFL': params[1],
        'sigma': xs
    }
