# SPDX-FileCopyrightText: 2026 OpenMC contributors
# SPDX-License-Identifier: MIT

from math import log

from pytest import approx

from endf.function import Tabulated1D


def test_integrate_linear_log():
    x0, x1 = 1., 4.
    y0, y1 = 2., 5.
    f = Tabulated1D([x0, x1], [y0, y1], [2], [3])

    slope = (y1 - y0)/log(x1/x0)

    def antiderivative(x):
        return y0*x + slope*(x*log(x/x0) - x)

    for a, b in ((x0, x1), (2., 3.)):
        expected = antiderivative(b) - antiderivative(a)
        assert f.integrate(a, b) == approx(expected)

    assert f.integral()[-1] == approx(antiderivative(x1) - antiderivative(x0))


def test_integrate_constant_log_linear():
    f = Tabulated1D([1., 4.], [2., 2.], [2], [4])

    assert f.integrate(1., 4.) == approx(6.)
    assert f.integrate(2., 3.) == approx(2.)
    assert f.integral()[-1] == approx(6.)


def test_integrate_inverse_log_log():
    f = Tabulated1D([1., 4.], [2., 0.5], [2], [5])

    assert f.integrate(1., 4.) == approx(2.*log(4.))
    assert f.integrate(2., 3.) == approx(2.*log(3./2.))
    assert f.integral()[-1] == approx(2.*log(4.))


def test_integrate_clipped_empty_range():
    f = Tabulated1D([1., 2.], [1., 1.])

    for a, b in ((-2., -1.), (-1., -2.), (3., 4.), (4., 3.)):
        assert f.integrate(a, b) == 0.


def test_integrate_clipped_partial_range():
    f = Tabulated1D([1., 2.], [1., 1.])

    assert f.integrate(0., 1.5) == approx(0.5)
    assert f.integrate(1.5, 0.) == approx(-0.5)
