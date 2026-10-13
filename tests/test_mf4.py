# SPDX-FileCopyrightText: 2026 OpenMC contributors
# SPDX-License-Identifier: MIT

import numpy as np
from numpy.polynomial.legendre import leggauss
import pytest
from endf.function import Tabulated1D
from endf.mf4 import AngleDistribution


@pytest.mark.parametrize('ltt', [1, 3])
def test_legendre_conversion(ltt):
    coefficients = [np.array([]), np.array([0.2, 0.1]), np.array([-0.1])]
    energy = np.array([1.e3, 1.e-5, 1.])
    tabulated = Tabulated1D([-1., 0., 1.], [0.25, 0.5, 0.75])
    data = {
        'LTT': ltt,
        'LI': 0,
        'legendre': {'E': energy, 'a_l': coefficients},
        'tabulated': {'E': np.array([1.e6]), 'mu': [tabulated]},
    }
    original_coefficients = [a_l.copy() for a_l in coefficients]
    distribution = AngleDistribution.from_dict(data)
    repeated = AngleDistribution.from_dict(data)

    expected_energy = [1.e3, 1.e-5, 1., 1.e6] if ltt == 3 else [1.e3, 1.e-5, 1.]
    np.testing.assert_array_equal(distribution.energy, expected_energy)
    np.testing.assert_array_equal(repeated.energy, expected_energy)
    np.testing.assert_array_equal(energy, [1.e3, 1.e-5, 1.])
    assert len(distribution.mu) == len(repeated.mu) == len(expected_energy)

    expected_coef = [[0.5], [0.5, 0.3, 0.25], [0.5, -0.15]]
    expected_values = [[0.5, 0.5, 0.5], [0.45, 0.375, 1.05], [0.65, 0.5, 0.35]]
    expected_moments = [[1., 0., 0.], [1., 0.2, 0.1], [1., -0.1, 0.]]
    cosine, weights = leggauss(8)
    basis = np.array([np.ones(8), cosine, (3 * cosine**2 - 1) / 2])
    for i, density in enumerate(distribution.mu[:3]):
        assert density.coef == pytest.approx(expected_coef[i], rel=0., abs=1.e-13)
        assert density([-1., 0., 1.]) == pytest.approx(expected_values[i], rel=0., abs=1.e-13)
        assert density(0.) == pytest.approx(expected_values[i][1], rel=0., abs=1.e-13)
        moments = basis @ (weights * density(cosine))
        assert moments == pytest.approx(expected_moments[i], rel=0., abs=1.e-13)
        np.testing.assert_array_equal(coefficients[i], original_coefficients[i])
        np.testing.assert_array_equal(density.coef, repeated.mu[i].coef)

    if ltt == 3:
        assert distribution.mu[-1] is repeated.mu[-1] is tabulated
        np.testing.assert_array_equal(data['tabulated']['E'], [1.e6])
        np.testing.assert_array_equal(tabulated.x, [-1., 0., 1.])
        np.testing.assert_array_equal(tabulated.y, [0.25, 0.5, 0.75])
        assert tabulated([-1., 0., 1.]) == pytest.approx([0.25, 0.5, 0.75])


def test_tabulated_conversion():
    energy = np.array([1.e6])
    mu = [Tabulated1D([-1., 1.], [0.4, 0.6])]
    data = {'LTT': 2, 'LI': 0, 'tabulated': {'E': energy, 'mu': mu}}

    distribution = AngleDistribution.from_dict(data)

    assert distribution.energy is energy
    assert distribution.mu is mu
