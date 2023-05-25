# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from abc import ABC
from typing import TextIO

import numpy as np

from .records import get_tab1_record, get_tab2_record, get_head_record


def parse_mf5(file_obj: TextIO) -> dict:
    ZA, AWR, _, _, NK, _ = get_head_record(file_obj)

    data = {'ZA': ZA, 'AWR': AWR, 'NK': NK}
    data['subsections'] = []
    for _ in range(NK):
        subsection = {}
        params, applicability = get_tab1_record(file_obj)
        subsection['LF'] = LF = params[3]
        subsection['p'] = applicability
        if LF == 1:
            dist = ArbitraryTabulated.dict_from_endf(file_obj, params)
        elif LF == 5:
            dist = GeneralEvaporation.dict_from_endf(file_obj, params)
        elif LF == 7:
            dist = MaxwellEnergy.dict_from_endf(file_obj, params)
        elif LF == 9:
            dist = Evaporation.dict_from_endf(file_obj, params)
        elif LF == 11:
            dist = WattEnergy.dict_from_endf(file_obj, params)
        elif LF == 12:
            dist = MadlandNix.dict_from_endf(file_obj, params)

        subsection['distribution'] = dist
        data['subsections'].append(subsection)

    return data



class EnergyDistribution(ABC):
    """Abstract superclass for all energy distributions."""
    def __init__(self):
        pass

    @staticmethod
    def from_endf(file_obj: TextIO, params: list):
        """Generate energy distribution from MF=5 data

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        A sub-class of :class:`EnergyDistribution`

        """
        LF = params[3]
        if LF == 1:
            return ArbitraryTabulated.from_endf(file_obj, params)
        elif LF == 5:
            return GeneralEvaporation.from_endf(file_obj, params)
        elif LF == 7:
            return MaxwellEnergy.from_endf(file_obj, params)
        elif LF == 9:
            return Evaporation.from_endf(file_obj, params)
        elif LF == 11:
            return WattEnergy.from_endf(file_obj, params)
        elif LF == 12:
            return MadlandNix.from_endf(file_obj, params)

    @staticmethod
    def from_dict(subsection: dict):
        LF = subsection['LF']
        data = subsection['distribution']
        if LF == 1:
            return ArbitraryTabulated.from_dict(data)
        elif LF == 5:
            return GeneralEvaporation.from_dict(data)
        elif LF == 7:
            return MaxwellEnergy.from_dict(data)
        elif LF == 9:
            return Evaporation.from_dict(data)
        elif LF == 11:
            return WattEnergy.from_dict(data)
        elif LF == 12:
            return MadlandNix.from_dict(data)


class ArbitraryTabulated(EnergyDistribution):
    r"""Arbitrary tabulated function given in ENDF MF=5, LF=1 represented as

    .. math::
         f(E \rightarrow E') = g(E \rightarrow E')

    Parameters
    ----------
    energy : numpy.ndarray
        Array of incident neutron energies
    pdf : list of openmc.data.Tabulated1D
        Tabulated outgoing energy distribution probability density functions

    Attributes
    ----------
    energy : numpy.ndarray
        Array of incident neutron energies
    pdf : list of openmc.data.Tabulated1D
        Tabulated outgoing energy distribution probability density functions

    """

    def __init__(self, energy, pdf):
        super().__init__()
        self.energy = energy
        self.pdf = pdf

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse arbitrary tabulated distribution (LF=1)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.

        Returns
        -------
        dict
            Arbitrary tabulated distribution data

        """
        data = {}
        params, data['E_int'] = get_tab2_record(file_obj)
        n_energies = params[5]

        energy = np.zeros(n_energies)
        pdf = []
        for j in range(n_energies):
            params, func = get_tab1_record(file_obj)
            energy[j] = params[1]
            pdf.append(func)
        data['E'] = energy
        data['g'] = pdf
        return data

    @classmethod
    def from_endf(cls, file_obj: TextIO, params: list):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['E'], data['g'])

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['E'], data['g'])



class GeneralEvaporation(EnergyDistribution):
    r"""General evaporation spectrum given in ENDF MF=5, LF=5 represented as

    .. math::
        f(E \rightarrow E') = g(E'/\theta(E))

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy :math:`E`
    g : openmc.data.Tabulated1D
        Tabulated function of :math:`x = E'/\theta(E)`
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy :math:`E`
    g : openmc.data.Tabulated1D
        Tabulated function of :math:`x = E'/\theta(E)`
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, g, u):
        super().__init__()
        self.theta = theta
        self.g = g
        self.u = u

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse general evaporation spectrum (MF=5)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        openmc.data.GeneralEvaporation
            General evaporation spectrum

        """
        _, theta = get_tab1_record(file_obj)
        _, g = get_tab1_record(file_obj)
        return {'U': params[0], 'theta': theta, 'g': g}

    @classmethod
    def from_endf(cls, file_obj: TextIO, params: list):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['theta'], data['g'], data['U'])

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['theta'], data['g'], data['U'])


class MaxwellEnergy(EnergyDistribution):
    r"""Simple Maxwellian fission spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{\sqrt{E'}}{I} e^{-E'/\theta(E)}

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        super().__init__()
        self.theta = theta
        self.u = u

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse Maxwellian fission spectrum (LF=7)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        dict
            Maxwellian distribution data

        """
        _, theta = get_tab1_record(file_obj)
        return {'U': params[0], 'theta': theta}

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['theta'], data['U'])


class Evaporation(EnergyDistribution):
    r"""Evaporation spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{E'}{I} e^{-E'/\theta(E)}

    Parameters
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    theta : openmc.data.Tabulated1D
        Tabulated function of incident neutron energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, theta, u):
        super().__init__()
        self.theta = theta
        self.u = u

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse evaporation spectrum (LF=9)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        data
            Evaporation spectrum data

        """
        _, theta = get_tab1_record(file_obj)
        return {'U': params[0], 'theta': theta}

    @classmethod
    def from_endf(cls, file_obj: TextIO, params: list):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['theta'], data['U'])

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['theta'], data['U'])


class WattEnergy(EnergyDistribution):
    r"""Energy-dependent Watt spectrum represented as

    .. math::
        f(E \rightarrow E') = \frac{e^{-E'/a}}{I} \sinh \left ( \sqrt{bE'}
        \right )

    Parameters
    ----------
    a, b : openmc.data.Tabulated1D
        Energy-dependent parameters tabulated as function of incident neutron
        energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    Attributes
    ----------
    a, b : openmc.data.Tabulated1D
        Energy-dependent parameters tabulated as function of incident neutron
        energy
    u : float
        Constant introduced to define the proper upper limit for the final
        particle energy such that :math:`0 \le E' \le E - U`

    """

    def __init__(self, a, b, u):
        super().__init__()
        self.a = a
        self.b = b
        self.u = u

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse energy-dependent Watt spectrum (MF=11)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        data
            Watt fission spectrum data

        """
        _, a = get_tab1_record(file_obj)
        _, b = get_tab1_record(file_obj)
        return {'U': params[0], 'a': a, 'b': b}

    @classmethod
    def from_endf(cls, file_obj: TextIO, params: list):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['a'], data['b'], data['U'])

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['a'], data['b'], data['U'])

class MadlandNix(EnergyDistribution):
    r"""Energy-dependent fission neutron spectrum (Madland and Nix) given in
    ENDF MF=5, LF=12 represented as

    .. math::
        f(E \rightarrow E') = \frac{1}{2} [ g(E', E_F(L)) + g(E', E_F(H))]

    where

    .. math::
        g(E',E_F) = \frac{1}{3\sqrt{E_F T_M}} \left [ u_2^{3/2} E_1 (u_2) -
        u_1^{3/2} E_1 (u_1) + \gamma \left ( \frac{3}{2}, u_2 \right ) - \gamma
        \left ( \frac{3}{2}, u_1 \right ) \right ] \\ u_1 = \left ( \sqrt{E'} -
        \sqrt{E_F} \right )^2 / T_M \\ u_2 = \left ( \sqrt{E'} + \sqrt{E_F}
        \right )^2 / T_M.

    Parameters
    ----------
    efl, efh : float
        Constants which represent the average kinetic energy per nucleon of the
        fission fragment (efl = light, efh = heavy)
    tm : openmc.data.Tabulated1D
        Parameter tabulated as a function of incident neutron energy

    Attributes
    ----------
    efl, efh : float
        Constants which represent the average kinetic energy per nucleon of the
        fission fragment (efl = light, efh = heavy)
    tm : openmc.data.Tabulated1D
        Parameter tabulated as a function of incident neutron energy

    """

    def __init__(self, efl, efh, tm):
        super().__init__()
        self.efl = efl
        self.efh = efh
        self.tm = tm

    @staticmethod
    def dict_from_endf(file_obj: TextIO, params: list) -> dict:
        """Parse Madland-Nix fission spectrum (LF=12)

        Parameters
        ----------
        file_obj : file-like object
            ENDF file positioned at the start of a section for an energy
            distribution.
        params : list
            List of parameters at the start of the energy distribution that
            includes the LF value indicating what type of energy distribution is
            present.

        Returns
        -------
        data
            Madland-Nix fission spectrum data

        """
        _, T_M = get_tab1_record(file_obj)
        return {'EFL': params[0], 'EFH': params[1], 'T_M': T_M}

    @classmethod
    def from_endf(cls, file_obj: TextIO, params: list):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['EFL'], data['EFH'], data['T_M'])

    @classmethod
    def from_dict(cls, data: dict):
        return cls(data['EFL'], data['EFH'], data['T_M'])


class LevelInelastic:
    r"""Level inelastic scattering

    Parameters
    ----------
    threshold : float
        Energy threshold in the laboratory system, :math:`(A + 1)/A * |Q|`
    mass_ratio : float
        :math:`(A/(A + 1))^2`

    Attributes
    ----------
    threshold : float
        Energy threshold in the laboratory system, :math:`(A + 1)/A * |Q|`
    mass_ratio : float
        :math:`(A/(A + 1))^2`

    """

    def __init__(self, threshold, mass_ratio):
        self.threshold = threshold
        self.mass_ratio = mass_ratio
