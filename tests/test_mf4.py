# SPDX-FileCopyrightText: 2026 OpenMC contributors
# SPDX-License-Identifier: MIT

import numpy as np
from numpy.polynomial import Legendre
from numpy.polynomial.legendre import leggauss
import pytest

from endf.function import Tabulated1D
from endf.mf4 import AngleDistribution


def angular_data(ltt, coefficients, energy=None):
    if energy is None:
        energy = np.arange(1, len(coefficients) + 1, dtype=float)
    data = {
        'LTT': ltt,
        'LI': 0,
        'legendre': {
            'E': np.asarray(energy),
            'a_l': [np.array(a_l, dtype=float) for a_l in coefficients],
        },
    }
    if ltt == 3:
        data['tabulated'] = {
            'E': np.array([1.e6, 2.e7]),
            'mu': [
                Tabulated1D([-1., 0., 1.], [0.25, 0.5, 0.75]),
                Tabulated1D([-1., 1.], [0.4, 0.6]),
            ],
        }
    return data


@pytest.mark.parametrize('ltt', [1, 3])
@pytest.mark.parametrize(
    'stored, expected_coef, expected_values',
    [
        ([], [0.5], [0.5, 0.5, 0.5]),
        ([0.2, 0.1], [0.5, 0.3, 0.25], [0.45, 0.375, 1.05]),
    ],
    ids=['isotropic', 'anisotropic'],
)
def test_legendre_values(ltt, stored, expected_coef, expected_values):
    distribution = AngleDistribution.from_dict(angular_data(ltt, [stored]))
    density = distribution.mu[0]

    assert isinstance(density, Legendre)
    np.testing.assert_allclose(density.coef, expected_coef, rtol=0., atol=1.e-13)
    np.testing.assert_allclose(
        density(np.array([-1., 0., 1.])), expected_values, rtol=0., atol=1.e-13
    )
    for cosine, expected in zip([-1., 0., 1.], expected_values):
        assert density(cosine) == pytest.approx(expected, rel=0., abs=1.e-13)


@pytest.mark.parametrize('ltt', [1, 3])
@pytest.mark.parametrize('stored', [[], [0.2, 0.1]], ids=['isotropic', 'anisotropic'])
def test_legendre_moments(ltt, stored):
    distribution = AngleDistribution.from_dict(angular_data(ltt, [stored]))
    cosine, weights = leggauss(8)
    density = distribution.mu[0](cosine)

    expected_moments = [1., *stored, *([0.] * (2 - len(stored)))]
    for order, expected in enumerate(expected_moments):
        moment = np.dot(weights, Legendre.basis(order)(cosine) * density)
        assert moment == pytest.approx(expected, rel=0., abs=1.e-13)


@pytest.mark.parametrize('ltt', [1, 3])
def test_legendre_multiple_energies_and_input_preservation(ltt):
    data = angular_data(ltt, [[], [0.2, 0.1], [-0.1]], [1.e3, 1.e-5, 1.])
    original_energy = data['legendre']['E'].copy()
    original_coefficients = [a_l.copy() for a_l in data['legendre']['a_l']]
    expected_energy = original_energy
    if ltt == 3:
        expected_energy = np.hstack((expected_energy, data['tabulated']['E']))

    first = AngleDistribution.from_dict(data)
    second = AngleDistribution.from_dict(data)
    np.testing.assert_array_equal(first.energy, expected_energy)
    np.testing.assert_array_equal(second.energy, expected_energy)
    np.testing.assert_array_equal(data['legendre']['E'], original_energy)
    assert len(first.mu) == len(second.mu) == len(expected_energy)

    cosine, weights = leggauss(8)
    for i, original in enumerate(original_coefficients):
        np.testing.assert_array_equal(data['legendre']['a_l'][i], original)
        np.testing.assert_array_equal(first.mu[i].coef, second.mu[i].coef)
        expected_moments = [1., *original, *([0.] * (2 - len(original)))]
        for order, expected in enumerate(expected_moments):
            moment = np.dot(
                weights, Legendre.basis(order)(cosine) * first.mu[i](cosine)
            )
            assert moment == pytest.approx(expected, rel=0., abs=1.e-13)


def test_mixed_tabulated_preservation():
    data = angular_data(3, [[], [0.2, 0.1]])
    tabulated = data['tabulated']
    original_energy = tabulated['E'].copy()
    original_values = [(f.x.copy(), f.y.copy()) for f in tabulated['mu']]

    distribution = AngleDistribution.from_dict(data)

    np.testing.assert_array_equal(distribution.energy[2:], original_energy)
    np.testing.assert_array_equal(tabulated['E'], original_energy)
    for actual, original, (x, y) in zip(
        distribution.mu[2:], tabulated['mu'], original_values
    ):
        assert actual is original
        np.testing.assert_array_equal(actual.x, x)
        np.testing.assert_array_equal(actual.y, y)
        np.testing.assert_array_equal(actual(x), y)


def test_pure_tabulated_preservation():
    tabulated = angular_data(3, [[]])['tabulated']
    data = {'LTT': 2, 'LI': 0, 'tabulated': tabulated}

    distribution = AngleDistribution.from_dict(data)

    assert distribution.energy is tabulated['E']
    assert distribution.mu is tabulated['mu']
