import numpy as np
from numpy.polynomial import Polynomial

from .function import Tabulated1D


class Product:
    """Secondary particle emitted in a nuclear reaction

    Parameters
    ----------
    name
        The particle type of the reaction product. Defaults to 'neutron'.
    yield_
        Yield of secondary particle in the reaction.
    distribution
        Distributions of energy and angle of product.
    applicability
        Probability of sampling a given distribution for this product.

    Attributes
    ----------
    applicability : Iterable of openmc.data.Tabulated1D
        Probability of sampling a given distribution for this product.
    decay_rate : float
        Decay rate in inverse seconds
    distribution : Iterable of AngleEnergy
        Distributions of energy and angle of product.
    emission_mode : {'prompt', 'delayed', 'total'}
        Indicate whether the particle is emitted immediately or whether it
        results from the decay of reaction product (e.g., neutron emitted from a
        delayed neutron precursor). A special value of 'total' is used when the
        yield represents particles from prompt and delayed sources.
    name : str
        The particle type of the reaction product
    yield_
        Yield of secondary particle in the reaction.

    """

    def __init__(self, name: str = 'neutron', yield_=None,
                 distribution=None, applicability=None):
        self.name = name
        if yield_ is None:
            self.yield_ = Polynomial((1,))  # 0-order polynomial, i.e., a constant
        else:
            self.yield_ = yield_
        self.decay_rate = 0.0
        if distribution is None:
            self.distribution = []
        else:
            self.distribution = distribution
        if applicability is None:
            self.applicability = []
        else:
            self.applicability = applicability
        self.emission_mode = 'prompt'

    def __repr__(self):
        if isinstance(self.yield_, Tabulated1D):
            if np.all(self.yield_.y == self.yield_.y[0]):
                return "<Product: {}, emission={}, yield={}>".format(
                    self.name, self.emission_mode, self.yield_.y[0])
            else:
                return "<Product: {}, emission={}, yield=tabulated>".format(
                    self.name, self.emission_mode)
        else:
            return "<Product: {}, emission={}, yield=polynomial>".format(
                self.name, self.emission_mode)
