# SPDX-FileCopyrightText: 2023-2025 Paul Romano
# SPDX-License-Identifier: MIT

from dataclasses import dataclass
from typing import Dict, List, Optional

from .function import Tabulated1D
from .material import Material

__all__ = ['RadionuclideProduction', 'radionuclide_production']


@dataclass
class RadionuclideProduction:
    """Production data for a single final state of one reaction.

    Joins the MF=8 identification of a radioactive product (its ZAP, LFS
    level number, and excitation energy) with the energy-dependent data
    evaluated for that state: the MF=9 yield multiplicity and/or the
    MF=10 production cross section. Note that LFS is a level index of
    the product nuclide, not an isomeric-state index; matching a level
    to a metastable state requires comparing :attr:`excitation_energy`
    against decay data.

    Attributes
    ----------
    ZAP : int
        1000*Z + A of the product nuclide
    LFS : int
        Level number of the final state (0 for the ground state)
    QM : float
        Mass-difference Q value in [eV]
    QI : float
        Reaction Q value for this state in [eV]
    ELFS : float or None
        Excitation energy of the final state in [eV] from MF=8, or None
        when the evaluation has no MF=8 subsection for this state
    yields : Tabulated1D or None
        MF=9 yield multiplicity of the reaction cross section as a
        function of incident energy in [eV]
    cross_section : Tabulated1D or None
        MF=10 production cross section in [b] as a function of incident
        energy in [eV]

    """

    ZAP: int
    LFS: int
    QM: float
    QI: float
    ELFS: Optional[float] = None
    yields: Optional[Tabulated1D] = None
    cross_section: Optional[Tabulated1D] = None

    @property
    def excitation_energy(self) -> float:
        """Excitation energy of the final state in [eV], taken from the
        MF=8 ELFS value when present and otherwise from QM minus QI."""
        if self.ELFS is not None:
            return self.ELFS
        return self.QM - self.QI


def radionuclide_production(
        material: Material) -> Dict[int, List[RadionuclideProduction]]:
    """Collect radionuclide production data from MF=8/9/10.

    For every reaction that has an MF=9 or MF=10 section, the final
    states are returned in the order they appear in the evaluation, with
    the MF=9 and MF=10 data for the same (ZAP, LFS) pair merged into one
    :class:`RadionuclideProduction` and the MF=8 excitation energy
    attached when available. The tabulated functions are returned
    exactly as evaluated.

    Parameters
    ----------
    material
        Material to read production data from

    Returns
    -------
    dict
        Mapping of MT numbers to lists of :class:`RadionuclideProduction`

    """
    by_mt: Dict[int, set] = {}
    for mf, mt in material.sections:
        if mf in (9, 10):
            by_mt.setdefault(mt, set()).add(mf)

    result = {}
    for mt in sorted(by_mt):
        # MF=8 links each (ZAP, LFS) pair to an excitation energy
        elfs = {}
        if (8, mt) in material.section_data:
            for subsection in material.section_data[8, mt]['subsections']:
                key = int(subsection['ZAP']), int(subsection['LFS'])
                elfs[key] = float(subsection['ELFS'])

        states: Dict[tuple, RadionuclideProduction] = {}
        ordered = []
        for mf in sorted(by_mt[mt]):
            for level in material.section_data[mf, mt]['levels']:
                key = int(level['IZAP']), int(level['LFS'])
                state = states.get(key)
                if state is None:
                    state = RadionuclideProduction(
                        ZAP=key[0],
                        LFS=key[1],
                        QM=level['QM'],
                        QI=level['QI'],
                        ELFS=elfs.get(key),
                    )
                    states[key] = state
                    ordered.append(state)
                if mf == 9:
                    state.yields = level['Y']
                else:
                    state.cross_section = level['sigma']
        result[mt] = ordered
    return result
