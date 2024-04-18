# SPDX-FileCopyrightText: Paul Romano
# SPDX-License-Identifier: MIT

from pytest import approx
from endf._records import float_endf


def test_float_sign():
    assert float_endf('+3.2146') == approx(3.2146)
    assert float_endf('-2.225002+6') == approx(-2.225002e6)


def test_float_no_leading_digit():
    assert float_endf('.12345') == approx(0.12345)


def test_float_double_digit_exponent():
    assert float_endf('6.022+23') == approx(6.022e23)
    assert float_endf('6.022-23') == approx(6.022e-23)


def test_float_whitespace():
    assert float_endf(' +1.01+ 2') == approx(101.0)
    assert float_endf(' -1.01- 2') == approx(-0.0101)
    assert float_endf('+ 2 . 3+ 1') == approx(23.0)
    assert float_endf('-7 .8 -1') == approx(-0.78)


def test_float_e_exponent():
    assert float_endf('3.14e0') == approx(3.14)
    assert float_endf('3.14E0') == approx(3.14)
    assert float_endf('3.14e-1') == approx(0.314)


def test_float_d_exponent():
    assert float_endf('3.14d0') == approx(3.14)
    assert float_endf('3.14D0') == approx(3.14)
    assert float_endf('3.14d-1') == approx(0.314)


def test_float_only_leading_digit():
    assert float_endf('1+2') == approx(100.0)
    assert float_endf('-1+2') == approx(-100.0)
    assert float_endf('1.+2') == approx(100.0)
    assert float_endf('-1.+2') == approx(-100.0)


def test_float_empty():
    assert float_endf('        ') == 0.0


def test_float_buffer_size():
    assert float_endf('9.876540000000000') == approx(9.87654)
