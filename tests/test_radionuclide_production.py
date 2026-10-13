from pathlib import Path

import pytest
import endf
from endf import RadionuclideProduction, radionuclide_production


@pytest.fixture
def in115():
    # ENDF/B-VIII.1 In115 evaluation trimmed to MF=1/451, the MF=3
    # sections for MT=4/16/102, and the corresponding MF=8/9/10 sections
    filename = Path(__file__).with_name('n-049_In-115_trimmed.endf')
    return endf.Material(filename)


def test_reactions_found(in115):
    production = radionuclide_production(in115)
    assert sorted(production) == [4, 16, 102]


def test_mf9_yields(in115):
    # In115 capture gives only the In116 first metastable state, with an
    # implicit ground-state share
    (state,) = radionuclide_production(in115)[102]
    assert state.ZAP == 49116
    assert state.LFS == 1
    assert state.ELFS == pytest.approx(127269.7)
    assert state.excitation_energy == pytest.approx(127269.7)
    assert state.cross_section is None
    assert state.yields is not None
    assert state.yields(0.0253) == pytest.approx(0.79)


def test_mf10_cross_section(in115):
    # Inelastic scattering to In115_m1 is given as an MF=10 partial
    # cross section with no MF=9 counterpart
    (state,) = radionuclide_production(in115)[4]
    assert state.ZAP == 49115
    assert state.LFS == 1
    assert state.ELFS == pytest.approx(336.2e3)
    assert state.yields is None
    assert state.cross_section is not None
    assert state.QM == 0.0
    assert state.QI == pytest.approx(-336.2e3)

    (state,) = radionuclide_production(in115)[16]
    assert state.ZAP == 49114
    assert state.LFS == 1
    assert state.ELFS == pytest.approx(190268.2)


def test_excitation_energy_fallback():
    # Without an MF=8 subsection the excitation energy falls back to the
    # difference of the Q values
    state = RadionuclideProduction(ZAP=41093, LFS=1, QM=0.0, QI=-30730.0)
    assert state.ELFS is None
    assert state.excitation_energy == pytest.approx(30730.0)


def test_material_without_data():
    filename = Path(__file__).with_name('n-095_Am_244.endf')
    material = endf.Material(filename)
    assert radionuclide_production(material) == {}
