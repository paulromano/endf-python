from typing import TextIO

from .records import get_head_record, get_tab1_record


def parse_mf3(file_obj: TextIO) -> dict:
    # Generate cross section
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
