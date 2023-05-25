# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from __future__ import annotations
from typing import TextIO
from warnings import warn

import numpy as np
from numpy.polynomial import Legendre

from .records import get_head_record, get_cont_record, get_tab2_record, \
    get_tab1_record, get_list_record


def parse_mf4(file_obj: TextIO) -> dict:
    # Read first two records
    ZA, AWR, LVT, LTT, _, _ = get_head_record(file_obj)
    _, _, LI, LCT, NK, NM = get_cont_record(file_obj)

    # initialize dictionary for angular distribution
    data = {'ZA': ZA, 'AWR': AWR, 'LTT': LTT, 'LI': LI, 'LCT': LCT}

    # Check for obsolete energy transformation matrix. If present, just skip
    # it and keep reading
    if LVT > 0:
        warn('Obsolete energy transformation matrix in MF=4 angular distribution.')
        for _ in range((NK + 5)//6):
            file_obj.readline()

    def legendre_data(file_obj):
        data = {}
        params, data['E_int'] = get_tab2_record(file_obj)
        n_energy = params[5]

        energy = np.zeros(n_energy)
        a_l = []
        for i in range(n_energy):
            items, al = get_list_record(file_obj)
            data['T'] = items[0]
            energy[i] = items[1]
            data['LT'] = items[2]
            coefficients = np.array(al)
            a_l.append(coefficients)
        data['a_l'] = a_l
        data['E'] = energy
        return data

    def tabulated_data(file_obj):
        data = {}
        params, data['E_int'] = get_tab2_record(file_obj)
        n_energy = params[5]

        energy = np.zeros(n_energy)
        mu = []
        for i in range(n_energy):
            params, f = get_tab1_record(file_obj)
            data['T'] = params[0]
            energy[i] = params[1]
            data['LT'] = params[2]
            mu.append(f)
        data['E'] = energy
        data['mu'] = mu
        return data

    if LTT == 0 and LI == 1:
        # Purely isotropic
        pass

    elif LTT == 1 and LI == 0:
        # Legendre polynomial coefficients
        data['legendre'] = legendre_data(file_obj)

    elif LTT == 2 and LI == 0:
        # Tabulated probability distribution
        data['tabulated'] = tabulated_data(file_obj)

    elif LTT == 3 and LI == 0:
        # Legendre for low energies / tabulated for high energies
        data['legendre'] = legendre_data(file_obj)
        data['tabulated'] = tabulated_data(file_obj)

    return data


class AngleDistribution:
    """Angle distribution as a function of incoming energy

    Parameters
    ----------
    energy
        Incoming energies in eV at which distributions exist
    mu
        Distribution of scattering cosines corresponding to each incoming energy

    Attributes
    ----------
    energy
        Incoming energies in eV at which distributions exist
    mu
        Distribution of scattering cosines corresponding to each incoming energy

    """

    def __init__(self, energy, mu):
        self.energy = energy
        self.mu = mu

    @classmethod
    def from_dict(cls, data: dict) -> AngleDistribution:
        LTT = data['LTT']
        LI = data['LI']
        if LTT == 0 and LI == 1:
            # Purely isotropic
            # TODO: Use uniform here
            energy = []
            mu = []
        elif LTT == 1 and LI == 0:
            energy = data['legendre']['E']
            mu = []
            for a_l in data['legendre']['a_l']:
                coef = np.insert(a_l, 0, 1.0)
                mu.append(Legendre(coef))
        elif LTT == 2 and LI == 0:
            energy = data['tabulated']['E']
            mu = data['tabulated']['mu']
        elif LTT == 3 and LI == 0:
            # Get Legendre first
            energy_leg = data['legendre']['E']
            mu_leg = []
            for a_l in data['legendre']['a_l']:
                coef = np.insert(a_l, 0, 1.0)
                mu_leg.append(Legendre(coef))

            # Then get tabulated
            energy_tab = data['tabulated']['E']
            mu_tab = data['tabulated']['mu']

            # Combine
            energy = np.hstack((energy_leg, energy_tab))
            mu = mu_leg + mu_tab

        return cls(energy, mu)

