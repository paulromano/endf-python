# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_head_record, get_list_record, get_tab1_record, get_cont_record


def parse_mf8(file_obj: TextIO) -> dict:
    """Parse radioactive nuclide production data from MF=8

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Radioactive nuclide production data
    """
    ZA, AWR, LIS, LISO, NS, NO = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LIS': LIS, 'LISO': LISO, 'NS': NS, 'NO': NO}

    data['subsections'] = []
    for _ in range(NS):
        if NO == 0:
            (ZAP, ELFS, LMF, LFS, _, _), values = get_list_record(file_obj)
            ND = len(values) // 6
            subsection = {'ZAP': ZAP, 'ELFS': ELFS, 'LMF': LMF, 'LFS': LFS, 'ND': ND}
            subsection['HL'] = values[::6]
            subsection['RTYP'] = values[1::6]
            subsection['ZAN'] = values[2::6]
            subsection['BR'] = values[3::6]
            subsection['END'] = values[4::6]
            subsection['CT'] = values[5::6]
        elif NO == 1:
            ZAP, ELFS, LMF, LFS, _, _ = get_cont_record(file_obj)
            subsection = {'ZAP': ZAP, 'ELFS': ELFS, 'LMF': LMF, 'LFS': LFS}
        data['subsections'].append(subsection)

    return data


def parse_mf8_mt454(file_obj: TextIO) -> dict:
    """Parse fission product yield data from MF=8, MT=454 / MT=459

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Fission product yield data

    """

    # Determine number of energies
    items = get_head_record(file_obj)
    data = {'ZA': items[0], 'AWR': items[1], 'LE': items[2] - 1}
    data['yields'] = []
    for i in range(data['LE'] + 1):
        # Determine i-th energy and number of products
        (E, _, I, _, NN, NFP), values = get_list_record(file_obj)
        yield_data = {'E': E, 'NN': NN, 'NFP': NFP}
        yield_data['LE' if i == 0 else 'I'] = I

        # Get yields for i-th energy
        yield_data['products'] = products = []
        for j in range(NFP):
            ZAFP = values[4*j]
            FPS = values[4*j + 1]
            Y = (values[4*j + 2], values[4*j + 3])
            products.append({'ZAFP': ZAFP, 'FPS': FPS, 'Y': Y})

        data['yields'].append(yield_data)

    return data


def parse_mf8_mt457(file_obj: TextIO) -> dict:
    """Parse radioactive decay data from MF=8, MT=457

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Radioactive decay data

    """
    # Get head record
    ZA, AWR, LIS, LISO, NST, NSP = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LIS': LIS, 'LISO': LISO, 'NST': NST, 'NSP': NSP}

    # Check if nuclide is stable
    if NST == 1:
        get_list_record(file_obj)
        (SPI, PAR, *_), values = get_list_record(file_obj)
        data['SPI'] = SPI
        data['PAR'] = PAR
        return data

    # Half-life and decay energies
    items, values = get_list_record(file_obj)
    data['T1/2'] = (items[0], items[1])
    data['NC'] = NC = items[4]//2
    data['Ex'] = list(zip(values[::2], values[1::2]))

    items, values = get_list_record(file_obj)
    data['SPI'], data['PAR'], *_ = items

    # Decay mode information
    data['NDK'] = NDK = items[5]  # Number of decay modes
    data['modes'] = []
    for i in range(NDK):
        RTYP = values[6*i]
        RFS = values[6*i + 1]
        Q = tuple(values[6*i + 2:6*i + 4])
        BR = tuple(values[6*i + 4:6*(i + 1)])
        mode = {'RTYP': RTYP, 'RFS': RFS, 'Q': Q, 'BR': BR}
        data['modes'].append(mode)

    # Read spectra
    data['spectra'] = []
    for i in range(NSP):
        items, values = get_list_record(file_obj)
        _, STYP, LCON, LCOV, _, NER = items
        spectrum = {'STYP': STYP, 'LCON': LCON, 'LCOV': LCOV, 'NER': NER}

        # Decay radiation type
        spectrum['FD'] = tuple(values[0:2])
        spectrum['ER_AV'] = tuple(values[2:4])
        spectrum['FC'] = tuple(values[4:6])

        if LCON != 1:
            # Information about discrete spectrum
            spectrum['discrete'] = []
            for j in range(NER):
                items, values = get_list_record(file_obj)
                discrete = {}
                discrete['ER'] = tuple(items[0:2])
                discrete['RTYP'] = values[0]
                discrete['TYPE'] = values[1]
                if STYP == 0:
                    discrete['RI'] = tuple(values[2:4])
                    discrete['RIS'] = tuple(values[4:6])
                    discrete['RICC'] = tuple(values[6:8])
                    discrete['RICK'] = tuple(values[8:10])
                    discrete['RICL'] = tuple(values[10:12])
                spectrum['discrete'].append(discrete)

        if LCON != 0:
            # Read continuous spectrum
            params, RP = get_tab1_record(file_obj)
            spectrum['continuous'] = {'RTYP': params[0], 'RP': RP}

        # Read continuous covariance (Ek, Fk) table
        if LCOV not in (0, 2) and LCON != 0:
            items, values = get_list_record(file_obj)
            covar_continuous = {'LB': items[3]}
            covar_continuous['Ek'] = np.array(values[::2])
            covar_continuous['Fk'] = np.array(values[1::2])
            spectrum['continuous_covariance'] = covar_continuous

        if LCOV not in (0, 1):
            (_, _, LS, LB, NE, NERP), values = get_list_record(file_obj)
            covar_discrete = {'LS': LS, 'LB': LB, 'NE': NE, 'NERP': NERP}
            covar_discrete['Ek'] = np.array(values[:NERP])
            covar_discrete['Fkk'] = np.array(values[NERP:])
            # TODO: Reorder and shape Fkk based on the packing order described
            # in section 8.4 of the ENDF manual
            spectrum['discrete_covariance'] = covar_discrete

        # Add spectrum to list
        data['spectra'].append(spectrum)

    return data
