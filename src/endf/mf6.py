# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_tab2_record, get_list_record, get_head_record, \
    get_tab1_record, get_cont_record


def parse_mf6(file_obj: TextIO) -> dict:
    """Parse product energy-angle distributions from MF=6

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Product energy-angle distribution data

    """
    # Read HEAD record
    ZA, AWR, JP, LCT, NK, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'JP': JP, 'LCT': LCT, 'NK': NK}

    data['products'] = products = []
    for i in range(NK):
        # Get yield for this product
        (ZAP, AWP, LIP, LAW), y_i = get_tab1_record(file_obj)
        ZAP = int(ZAP)

        p = {'ZAP': ZAP, 'AWP': AWP, 'LIP': LIP, 'LAW': LAW, 'y_i': y_i}

        if LAW < 0:
            # Distribution given elsewhere
            pass
        elif LAW == 0:
            # No distribution given
            pass
        elif LAW == 1:
            # Continuum energy-angle distribution
            p['distribution'] = ContinuumEnergyAngle.dict_from_endf(file_obj)

        elif LAW == 2:
            # Discrete two-body scattering
            p['distribution'] = DiscreteTwoBodyScattering.dict_from_endf(file_obj)
        elif LAW == 3:
            # Isotropic discrete emission
            pass

        elif LAW == 4:
            # Discrete two-body recoil
            pass

        elif LAW == 5:
            # Charged particle elastic scattering
            p['distribution'] = ChargedParticleElasticScattering.dict_from_endf(file_obj)

        elif LAW == 6:
            # N-body phase-space distribution
            p['distribution'] = NBodyPhaseSpace.dict_from_endf(file_obj)

        elif LAW == 7:
            # Laboratory energy-angle distribution
            p['distribution'] = LaboratoryAngleEnergy.dict_from_endf(file_obj)

        products.append(p)

    return data


class ContinuumEnergyAngle:
    def __init__(self):
        pass

    @staticmethod
    def dict_from_endf(file_obj: TextIO) -> dict:
        params, E_int = get_tab2_record(file_obj)
        _, _, LANG, LEP, NR, NE = params

        data = {'LANG': LANG, 'LEP': LEP, 'NR': NR, 'NE': NE, 'E_int': E_int}

        data['E'] = np.zeros(NE)
        data['distribution'] = []
        for i in range(NE):
            items, values = get_list_record(file_obj)
            _, E_i, ND, NA, NW, NEP = items
            dist = {'ND': ND, 'NA': NA, 'NW': NW, 'NEP': NEP}
            data['E'][i] = E_i
            values = np.asarray(values)
            values.shape = (NEP, NA + 2)
            dist["E'"] = values[:, 0]
            dist['b'] = values[:, 1:]
            data['distribution'].append(dist)

        return data


class UncorrelatedAngleEnergy:
    def __init__(self):
        pass


class DiscreteTwoBodyScattering:
    def __init__(self):
        pass

    @staticmethod
    def dict_from_endf(file_obj: TextIO) -> dict:
        params, E_int = get_tab2_record(file_obj)
        *_, NR, NE = params
        data = {'NR': NR, 'NE': NE, 'E_int': E_int}

        data['E'] = np.zeros(NE)
        data['distribution'] = []
        for i in range(NE):
            items, values = get_list_record(file_obj)
            _, E_i, LANG, _, NW, NL = items
            dist = {'LANG': LANG, 'NW': NW, 'NL': NL}
            data['E'][i] = E_i
            data['A_l'] = np.asarray(values)
            data['distribution'].append(dist)


class ChargedParticleElasticScattering:
    def __init__(self):
        pass

    @staticmethod
    def dict_from_endf(file_obj: TextIO) -> dict:
        (SPI, _, LIDP, _, NR, NE), E_int = get_tab2_record(file_obj)
        data = {'SPI': SPI, 'LIDP': LIDP, 'NE': NE, 'E_int': E_int}

        # Read distribution data for each incident energy
        data['distribution'] = []
        for _ in range(NE):
            (_, E, LTP, _, NW, NL), A = get_list_record(file_obj)
            dist = {'E': E, 'LTP': LTP, 'NW': NW, 'NL': NL, 'A': A}
            data['distribution'].append(dist)

        return data


class NBodyPhaseSpace:
    def __init__(self):
        pass

    @staticmethod
    def dict_from_endf(file_obj: TextIO) -> dict:
        APSX, *_, NPSX = get_cont_record(file_obj)
        return {'APSX': APSX, 'NPSX': NPSX}


class LaboratoryAngleEnergy:
    def __init__(self):
        pass

    @staticmethod
    def dict_from_endf(file_obj: TextIO) -> dict:
        # Read top-level TAB2 record
        (*_, NR, NE), E_int = get_tab2_record(file_obj)
        data = {'NE': NE, 'E_int': E_int}

        data['distribution'] = []
        for _ in range(NE):
            # Read TAB2 record for the i-th incident energy
            (_, E, _, _, NRM, NMU), mu_int = get_tab2_record(file_obj)
            dist = {'E': E, 'NRM': NRM, 'NMU': NMU, 'mu_int': mu_int}

            dist['mu'] = []
            for _ in range(NMU):
                # Read TAB1 record for the j-th outgoing cosine
                (_, mu, *_), f = get_tab1_record(file_obj)
                dist['mu'].append({'mu': mu, 'f': f})

            data['distribution'].append(dist)

        return data
