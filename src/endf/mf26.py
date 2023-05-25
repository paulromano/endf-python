# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO
from warnings import warn

from .mf6 import ContinuumEnergyAngle, DiscreteTwoBodyScattering
from .records import get_head_record, get_tab1_record


def parse_mf26(file_obj: TextIO) -> dict:
    """Parse secondary distributions for atomic data from MF=26

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Secondary distribution data

    """
    ZA, AWR, _, _, NK, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'NK': NK}

    data['products'] = products = []
    for i in range(NK):
        # Get yield for this product
        (ZAP, AWI, _, LAW), y = get_tab1_record(file_obj)
        ZAP = int(ZAP)

        p = {'ZAP': ZAP, 'AWI': AWI, 'LAW': LAW, 'y': y}

        if LAW == 1:
            # Continuum energy-angle distribution
            p['distribution'] = ContinuumEnergyAngle.dict_from_endf(file_obj)

        elif LAW == 2:
            # Discrete two-body scattering
            p['distribution'] = DiscreteTwoBodyScattering.dict_from_endf(file_obj)

        elif LAW == 8:
            # Energy transfer for excitation
            _, ET = get_tab1_record(file_obj)
            p['distribution'] = {'ET': ET}

        else:
            warn(f'Unrecognized {LAW=} in MF=26')

        products.append(p)

    return data
