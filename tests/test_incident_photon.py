# SPDX-FileCopyrightText: Paul Romano
# SPDX-License-Identifier: MIT

from pathlib import Path

import numpy as np
import pytest
import endf
from endf.incident_photon import (
    AtomicRelaxation, PhotonReaction, PHOTON_REACTION_NAME, PHOTON_REACTION_MT
)


@pytest.fixture
def h_photoatomic():
    filename = Path(__file__).with_name('photoat-001_H_000.endf')
    return endf.Material(filename)


@pytest.fixture
def h_relaxation():
    filename = Path(__file__).with_name('atom-001_H_000.endf')
    return endf.Material(filename)


@pytest.fixture
def hydrogen(h_photoatomic, h_relaxation):
    return endf.IncidentPhoton.from_endf(h_photoatomic, relaxation=h_relaxation)


def test_from_endf(hydrogen):
    assert isinstance(hydrogen, endf.IncidentPhoton)
    assert hydrogen.atomic_number == 1
    assert hydrogen.name == 'H'


def test_reactions(hydrogen):
    assert len(hydrogen.reactions) == 8
    expected_mts = {501, 502, 504, 515, 516, 517, 522, 534}
    assert set(hydrogen.reactions.keys()) == expected_mts
    for rx in hydrogen:
        assert isinstance(rx, PhotonReaction)
        assert rx.xs is not None


def test_contains(hydrogen):
    assert 502 in hydrogen
    assert 504 in hydrogen
    assert 999 not in hydrogen


def test_getitem_by_mt(hydrogen):
    rx = hydrogen[502]
    assert isinstance(rx, PhotonReaction)
    assert rx.MT == 502


def test_getitem_by_name(hydrogen):
    assert hydrogen['coherent'].MT == 502
    assert hydrogen['incoherent'].MT == 504
    assert hydrogen['photoelectric'].MT == 522
    assert hydrogen['total'].MT == 501
    assert hydrogen['K'].MT == 534


def test_getitem_invalid(hydrogen):
    with pytest.raises(ValueError, match="No reaction with label"):
        hydrogen['nonexistent']
    with pytest.raises(ValueError, match="No reaction with"):
        hydrogen[999]


def test_repr(hydrogen):
    r = repr(hydrogen)
    assert 'IncidentPhoton' in r
    assert 'H' in r
    assert '8 reactions' in r


def test_iter(hydrogen):
    reactions = list(hydrogen)
    assert len(reactions) == 8
    assert all(isinstance(rx, PhotonReaction) for rx in reactions)


def test_coherent_scattering_factor(hydrogen):
    rx = hydrogen[502]
    assert rx.scattering_factor is not None
    assert rx.scattering_factor.x.size == 1253


def test_anomalous_scattering(hydrogen):
    rx = hydrogen[502]
    assert rx.anomalous_real is not None
    assert rx.anomalous_imag is not None
    assert rx.anomalous_real.x.size == 297
    assert rx.anomalous_imag.x.size == 297


def test_incoherent_scattering_factor(hydrogen):
    rx = hydrogen[504]
    assert rx.scattering_factor is not None
    assert rx.scattering_factor.x.size == 398


def test_subshell_reaction(hydrogen):
    rx = hydrogen[534]
    assert rx.subshell_binding_energy == pytest.approx(13.6)
    assert rx.fluorescence_yield == pytest.approx(0.0)


def test_atomic_relaxation(hydrogen):
    ar = hydrogen.atomic_relaxation
    assert isinstance(ar, AtomicRelaxation)
    assert ar.subshells == ['K']
    assert ar.binding_energy['K'] == pytest.approx(13.6)
    assert ar.num_electrons['K'] == pytest.approx(1.0)


def test_atomic_relaxation_repr(hydrogen):
    assert '1 subshells' in repr(hydrogen.atomic_relaxation)


def test_from_endf_without_relaxation(h_photoatomic):
    ip = endf.IncidentPhoton.from_endf(h_photoatomic)
    assert ip.atomic_number == 1
    assert ip.atomic_relaxation is None
    assert 502 in ip


def test_from_endf_file_path():
    filename = Path(__file__).with_name('photoat-001_H_000.endf')
    ip = endf.IncidentPhoton.from_endf(filename)
    assert ip.atomic_number == 1
    assert len(ip.reactions) == 8


def test_interpret(h_photoatomic):
    ip = h_photoatomic.interpret()
    assert isinstance(ip, endf.IncidentPhoton)
    assert ip.atomic_number == 1


def test_photon_reaction_repr():
    rx = PhotonReaction(502)
    assert 'MT=502' in repr(rx)
    assert 'coherent' in repr(rx)

    rx2 = PhotonReaction(9999)
    assert 'MT=9999' in repr(rx2)


def test_reaction_name_mapping():
    assert PHOTON_REACTION_NAME[502] == 'coherent'
    assert PHOTON_REACTION_NAME[504] == 'incoherent'
    assert PHOTON_REACTION_NAME[534] == 'K'
    assert PHOTON_REACTION_MT['coherent'] == 502
    assert PHOTON_REACTION_MT['K'] == 534


def test_get_reaction_components(hydrogen):
    # MT=522 (photoelectric) should decompose to subshell reactions
    components = hydrogen._get_reaction_components(522)
    assert 534 in components

    # MT=501 (total) should decompose
    components = hydrogen._get_reaction_components(501)
    assert len(components) > 1

    # A non-redundant reaction returns itself
    components = hydrogen._get_reaction_components(502)
    assert components == [502]

    # A missing reaction returns empty
    components = hydrogen._get_reaction_components(999)
    assert components == []
