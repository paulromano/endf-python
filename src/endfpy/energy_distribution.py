from abc import ABC

import numpy as np

from .records import get_tab1_record, get_tab2_record


class EnergyDistribution(ABC):
    """Abstract superclass for all energy distributions."""
    def __init__(self):
        pass

    @staticmethod
    def from_endf(file_obj, params):
        """Generate energy distribution from an ENDF evaluation

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
        openmc.data.EnergyDistribution
            A sub-class of :class:`openmc.data.EnergyDistribution`

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
    def dict_from_endf(file_obj, params):
        """Generate arbitrary tabulated distribution from an ENDF evaluation

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
        openmc.data.ArbitraryTabulated
            Arbitrary tabulated distribution

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
    def from_endf(cls, file_obj, params):
        data = cls.dict_from_endf(file_obj, params)
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
    def dict_from_endf(file_obj, params):
        """Generate general evaporation spectrum from an ENDF evaluation

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
    def from_endf(cls, file_obj, params):
        data = cls.dict_from_endf(file_obj, params)
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
    def dict_from_endf(file_obj, params):
        """Generate Maxwell distribution from an ENDF evaluation

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
        openmc.data.MaxwellEnergy
            Maxwell distribution

        """
        _, theta = get_tab1_record(file_obj)
        return {'U': params[0], 'theta': theta}

    @classmethod
    def from_endf(cls, file_obj, params):
        data = cls.dict_from_endf(file_obj, params)
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
    def dict_from_endf(file_obj, params):
        """Generate evaporation spectrum from an ENDF evaluation

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
        openmc.data.Evaporation
            Evaporation spectrum

        """
        _, theta = get_tab1_record(file_obj)
        return {'U': params[0], 'theta': theta}

    @classmethod
    def from_endf(cls, file_obj, params):
        data = cls.dict_from_endf(file_obj, params)
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
    def dict_from_endf(file_obj, params):
        """Generate Watt fission spectrum from an ENDF evaluation

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
        openmc.data.WattEnergy
            Watt fission spectrum

        """
        _, a = get_tab1_record(file_obj)
        _, b = get_tab1_record(file_obj)
        return {'U': params[0], 'a': a, 'b': b}


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
    def dict_from_endf(file_obj, params):
        """Generate Madland-Nix fission spectrum from an ENDF evaluation

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
        openmc.data.MadlandNix
            Madland-Nix fission spectrum

        """
        _, T_M = get_tab1_record(file_obj)
        return {'EFL': params[0], 'EFH': params[1], 'T_M': T_M}

    @classmethod
    def from_endf(cls, file_obj, params):
        data = cls.dict_from_endf(file_obj, params)
        return cls(data['EFL'], data['EFH'], data['T_M'])
