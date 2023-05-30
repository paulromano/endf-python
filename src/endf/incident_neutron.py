# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from __future__ import annotations
from typing import Union, List

import numpy as np

from .data import gnds_name, temperature_str, ATOMIC_SYMBOL, EV_PER_MEV
from .material import Material
from .fileutils import PathLike
from .function import Tabulated1D
from .reaction import Reaction, REACTION_MT
from . import ace


SUM_RULES = {
    1: [2, 3],
    3: [4, 5, 11, 16, 17, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35,
        36, 37, 41, 42, 44, 45, 152, 153, 154, 156, 157, 158, 159, 160,
        161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
        173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185,
        186, 187, 188, 189, 190, 194, 195, 196, 198, 199, 200],
    4: list(range(50, 92)),
    16: list(range(875, 892)),
    18: [19, 20, 21, 38],
    27: [18, 101],
    101: [102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114,
        115, 116, 117, 155, 182, 191, 192, 193, 197],
    103: list(range(600, 650)),
    104: list(range(650, 700)),
    105: list(range(700, 750)),
    106: list(range(750, 800)),
    107: list(range(800, 850))
}


class IncidentNeutron:
    """Continuous-energy neutron interaction data.

    This class stores data derived from an ENDF-6 format neutron interaction
    sublibrary.

    Parameters
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates the
        ground state.

    Attributes
    ----------
    atomic_number : int
        Number of protons in the target nucleus
    atomic_symbol : str
        Atomic symbol of the nuclide, e.g., 'Zr'
    mass_number : int
        Number of nucleons in the target nucleus
    metastable : int
        Metastable state of the target nucleus. A value of zero indicates the
        ground state.
    name : str
        Name of the nuclide using the GNDS naming convention
    reactions : dict
        Contains the cross sections, secondary angle and energy distributions,
        and other associated data for each reaction. The keys are the MT values
        and the values are Reaction objects.
    """

    def __init__(self, atomic_number: int, mass_number: int, metastable: int = 0):
        self.atomic_number = atomic_number
        self.mass_number = mass_number
        self.metastable = metastable
        self.reactions = {}

    @classmethod
    def from_endf(cls, filename_or_mat: Union[PathLike, Material]) -> IncidentNeutron:
        """Generate incident neutron data from an ENDF file

        Parameters
        ----------
        filename_or_mat
            Path to ENDF-6 formatted file or material object

        Returns
        -------
        Incident neutron data

        """
        if not isinstance(filename_or_mat, Material):
            material = Material(filename_or_mat)
        else:
            material = filename_or_mat

        # Determine atomic number, mass number, and metastable state
        metadata = material[1, 451]
        Z, A = divmod(metadata['ZA'], 1000)
        data = cls(Z, A, metadata['LISO'])

        # Read each reaction
        for MF, MT in material.sections:
            if MF == 3:
                data.reactions[MT] = Reaction.from_endf(MT, material)
        return data

    @classmethod
    def from_ace(
        cls,
        filename_or_table: Union[PathLike, ace.Table],
        metastable_scheme: str = 'mcnp'
    ) -> IncidentNeutron:
        """Generate incident neutron continuous-energy data from an ACE table

        Parameters
        ----------
        ace_or_filename
            ACE table to read from. If the value is a string, it is assumed to
            be the filename for the ACE file.
        metastable_scheme : {'mcnp', 'nndc'}
            Determine how ZAID identifiers are to be interpreted in the case of
            a metastable nuclide. Because the normal ZAID (=1000*Z + A) does not
            encode metastable information, different conventions are used among
            different libraries. In MCNP libraries, the convention is to add 400
            for a metastable nuclide except for Am242m, for which 95242 is
            metastable and 95642 (or 1095242 in newer libraries) is the ground
            state. For NNDC libraries, ZAID is given as 1000*Z + A + 100*m.

        Returns
        -------
        Incident neutron continuous-energy data
        """
        # First obtain the data for the first provided ACE table/file
        if isinstance(filename_or_table, ace.Table):
            table = filename_or_table
        else:
            table = ace.get_table(filename_or_table)

        # If mass number hasn't been specified, make an educated guess
        zaid, xs = table.name.split('.')
        if not xs.endswith('c'):
            raise TypeError(f"{table} is not a continuous-energy neutron ACE table.")
        name, _, Z, mass_number, metastable = \
            ace.get_metadata(int(zaid), metastable_scheme)

        # Get string of temperature to use as a dictionary key
        strT = temperature_str(table.temperature)

        # Create IncidentNeutron object (reactions will be added after)
        data = cls(Z, mass_number, metastable)

        # Read energy grid
        n_energy = table.nxs[3]
        i = table.jxs[1]
        energy = table.xss[i : i + n_energy]*EV_PER_MEV
        total_xs = table.xss[i + n_energy : i + 2*n_energy]
        absorption_xs = table.xss[i + 2*n_energy : i + 3*n_energy]
        heating_number = table.xss[i + 4*n_energy : i + 5*n_energy]*EV_PER_MEV

        # Create redundant reaction for total (MT=1)
        xs = {strT: Tabulated1D(energy, total_xs)}
        data.reactions[1] = Reaction(1, xs, redundant=True)

        # Create redundant reaction for absorption (MT=101)
        if np.count_nonzero(absorption_xs) > 0:
            xs = {strT: Tabulated1D(energy, absorption_xs)}
            data.reactions[101] = Reaction(101, xs, redundant=True)

        # Create redundant reaction for heating (MT=301)
        xs = {strT: Tabulated1D(energy, heating_number*total_xs)}
        data.reactions[301] = Reaction(301, xs, redundant=True)

        # Read each reaction
        n_reaction = table.nxs[4] + 1
        for i in range(n_reaction):
            rx = Reaction.from_ace(table, i)
            data.reactions[rx.MT] = rx

        # Make sure redundant cross sections that are present in an ACE file get
        # marked as such
        for rx in data:
            mts = data._get_reaction_components(rx.MT)
            if mts != [rx.MT]:
                rx.redundant = True
            if rx.MT in (203, 204, 205, 206, 207, 444):
                rx.redundant = True

        return data

    def __contains__(self, MT: int):
        return MT in self.reactions

    def __getitem__(self, MT_or_name: int) -> Reaction:
        if isinstance(MT_or_name, str):
            if MT_or_name in REACTION_MT:
                MT = REACTION_MT[MT_or_name]
            elif f'({MT_or_name})' in REACTION_MT:
                MT = REACTION_MT[f'({MT_or_name})']
            else:
                raise ValueError(f"No reaction with label {MT_or_name}")
        else:
            MT = MT_or_name

        if MT in self.reactions:
            return self.reactions[MT]
        else:
            # TODO: Try to create a redundant cross section
            raise ValueError(f"No reaction with {MT=}")

    def __repr__(self) -> str:
        return f"<IncidentNeutron: {self.name}, {len(self.reactions)} reactions>"

    def __iter__(self):
        return iter(self.reactions.values())

    @property
    def name(self) -> str:
        return gnds_name(self.atomic_number, self.mass_number, self.metastable)

    @property
    def atomic_symbol(self) -> str:
        return ATOMIC_SYMBOL[self.atomic_number]

    def _get_reaction_components(self, MT: int) -> List[int]:
        """Determine what reactions make up redundant reaction.

        Parameters
        ----------
        mt : int
            ENDF MT number of the reaction to find components of.

        Returns
        -------
        mts : list of int
            ENDF MT numbers of reactions that make up the redundant reaction and
            have cross sections provided.

        """
        mts = []
        if MT in SUM_RULES:
            for MT_i in SUM_RULES[MT]:
                mts += self._get_reaction_components(MT_i)
        if mts:
            return mts
        else:
            return [MT] if MT in self else []
