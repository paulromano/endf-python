# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_list_record


def parse_mf28(file_obj: TextIO) -> dict:
    """Parse atomic relaxation data from MF=27

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Atomic relaxation data

    """
    ZA, AWR, _, _, NSS, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'NSS': NSS, 'shells': []}
    for _ in range(NSS):
        # Read LIST record for each shell
        (SUBI, _, _, _, NW, NTR), values = get_list_record(file_obj)
        shell = {'SUBI': SUBI, 'NTR': NTR}
        shell['EBI'] = values[0]
        shell['ELN'] = values[1]
        shell['SUBJ'] = values[6::6]
        shell['SUBK'] = values[7::6]
        shell['ETR'] = values[8::6]
        shell['FTR'] = values[9::6]
        data['shells'].append(shell)
    return data
