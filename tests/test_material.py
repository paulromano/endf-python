# SPDX-FileCopyrightText: Paul Romano
# SPDX-License-Identifier: MIT

from pathlib import Path

import pytest
import endf


@pytest.fixture
def am244():
    filename = Path(__file__).with_name('n-095_Am_244.endf')
    return endf.Material(filename)


def test_init_file():
    filename = Path(__file__).with_name('n-095_Am_244.endf')
    with open(filename) as fh:
        am244 = endf.Material(fh)
    assert am244.MAT == 9552


def test_sections(am244):
    assert isinstance(am244.sections, list)
    for mf_mt in am244.sections:
        assert len(mf_mt) == 2


def test_section_text(am244):
    for key, value in am244.section_text.items():
        assert isinstance(key, tuple)
        assert isinstance(value, str)


def test_section_data(am244):
    # Spot check metadata in MF=1, MT=451
    metadata = am244.section_data[1, 451]
    assert metadata['ZSYMAM'].strip() == '95-Am-244'
    assert metadata['EMAX'] == pytest.approx(20.0e6)

    # Spot check cross section data in MF=4
    capture = am244.section_data[3, 102]
    assert capture['QM'] == pytest.approx(6052990.0)

    # Indexing Material directly is same as indexing section_data
    assert am244[3, 102] is am244.section_data[3, 102]


def test_contains(am244):
    assert (3, 102) in am244
    assert (102, 3) not in am244


def test_repr(am244):
    assert '95-Am-244' in repr(am244)
    assert 'ENDF/B' in repr(am244)


def test_interpret(am244):
    am244_high_level = am244.interpret()
    assert isinstance(am244_high_level, endf.IncidentNeutron)
