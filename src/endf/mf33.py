# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

from .records import get_head_record, get_cont_record, get_list_record


def parse_mf33_subsection(file_obj) -> dict:
    XMF1, XLFS1, MAT1, MT1, NC, NI = get_cont_record(file_obj)
    subsection = {'XMF1': XMF1, 'XLFS1': XLFS1, 'MAT1': MAT1, 'MT1': MT1,
                  'NC': NC, 'NI': NI}

    subsection['nc_subsections'] = []
    for _ in range(NC):
        LTY = get_cont_record(file_obj)[3]
        if LTY == 0:
            (E1, E2, *_, NCI), values = get_list_record(file_obj)
            subsub = {'LTY': LTY, 'E1': E1, 'E2': E2, 'NCI': NCI}
            subsub['CI'] = values[::2]
            subsub['XMTI'] = values[1::2]
            subsection['nc_subsections'].append(subsub)
        else:
            (E1, E2, MATS, MTS, _, NEI), values = get_list_record(file_obj)
            subsub = {'LTY': LTY, 'E1': E1, 'E2': E2, 'MATS': MATS,
                        'MTS': MTS, 'NEI': NEI}
            subsub['XMFS'] = values[0]
            subsub['XLFSS'] = values[1]
            subsub['EI'] = values[2::2]
            subsub['WEI'] = values[3::2]
        subsection['nc_subsections'].append(subsub)

    subsection['ni_subsections'] = []
    for _ in range(NI):
        # Look ahead to determine LB
        pos = file_obj.tell()
        LB = get_cont_record(file_obj)[3]
        file_obj.seek(pos)

        if 0 <= LB <= 4:
            (_, _, LT, LB, NT, NP), values = get_list_record(file_obj)
            subsub = {'LT': LT, 'LB': LB, 'NT': NT, 'NP': NP}
            k_array = values[:NT - NP]
            subsub['Ek'] = k_array[::2]
            subsub['Fk'] = k_array[1::2]
            l_array = values[NT - NP:]
            subsub['El'] = l_array[::2]
            subsub['Fl'] = l_array[1::2]
        elif LB == 5:
            (_, _, LS, LB, NT, NE), values = get_list_record(file_obj)
            subsub = {'LS': LS, 'LB': LB, 'NT': NT, 'NE': NE}
            subsub['Ek'] = values[:NE]
            # TODO: Reoder/reshape values for Fk,k' matrix
            subsub['Fkk'] = values[NE:]
        elif LB == 6:
            (_, _, _, LB, NT, NER), values = get_list_record(file_obj)
            NEC = (NT - 1)//NER
            subsub = {'LB': LB, 'NT': NT, 'NER': NER, 'NEC': NEC}
            subsub['ER'] = values[:NER]
            subsub['EC'] = values[NER:NER + NEC]
            # TODO: Reorder/reshape values for Fkl matrix
            subsub['Fkl'] = values[NER + NEC:]
        elif LB in (8, 9):
            (_, _, LT, LB, NT, NP), values = get_list_record(file_obj)
            subsub = {'LT': LT, 'LB': LB, 'NT': NT, 'NP': NP}
            subsub['Ek'] = values[::2]
            subsub['Fk'] = values[1::2]
        else:
            raise ValueError(f"Unrecognized {LB=}")

        subsection['ni_subsections'].append(subsub)

    return subsection


def parse_mf33(file_obj: TextIO) -> dict:
    """Parse covariances of neutron cross sections from MF=33

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Cross section covariance data

    """
    ZA, AWR, _, MTL, _, NL = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'MTL': MTL, 'NL': NL, 'subsections': []}
    for _ in range(NL):
        subsection = parse_mf33_subsection(file_obj)
        data['subsections'].append(subsection)

    return data
