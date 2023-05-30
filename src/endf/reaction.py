# SPDX-FileCopyrightText: 2023 OpenMC contributors and Paul Romano
# SPDX-License-Identifier: MIT

from __future__ import annotations
from copy import deepcopy
from typing import List, Tuple
from warnings import warn

import numpy as np
from numpy.polynomial import Polynomial

from .data import gnds_name, temperature_str, ATOMIC_SYMBOL, EV_PER_MEV
from .material import Material
from .function import Tabulated1D
from .mf4 import AngleDistribution
from .mf5 import EnergyDistribution, LevelInelastic
from .mf6 import UncorrelatedAngleEnergy
from .product import Product
from . import ace


REACTION_NAME = {
    1: '(n,total)', 2: '(n,elastic)', 4: '(n,level)',
    5: '(n,misc)', 11: '(n,2nd)', 16: '(n,2n)', 17: '(n,3n)',
    18: '(n,fission)', 19: '(n,f)', 20: '(n,nf)', 21: '(n,2nf)',
    22: '(n,na)', 23: '(n,n3a)', 24: '(n,2na)', 25: '(n,3na)',
    27: '(n,absorption)', 28: '(n,np)', 29: '(n,n2a)',
    30: '(n,2n2a)', 32: '(n,nd)', 33: '(n,nt)', 34: '(n,n3He)',
    35: '(n,nd2a)', 36: '(n,nt2a)', 37: '(n,4n)', 38: '(n,3nf)',
    41: '(n,2np)', 42: '(n,3np)', 44: '(n,n2p)', 45: '(n,npa)',
    91: '(n,nc)', 101: '(n,disappear)', 102: '(n,gamma)',
    103: '(n,p)', 104: '(n,d)', 105: '(n,t)', 106: '(n,3He)',
    107: '(n,a)', 108: '(n,2a)', 109: '(n,3a)', 111: '(n,2p)',
    112: '(n,pa)', 113: '(n,t2a)', 114: '(n,d2a)', 115: '(n,pd)',
    116: '(n,pt)', 117: '(n,da)', 152: '(n,5n)', 153: '(n,6n)',
    154: '(n,2nt)', 155: '(n,ta)', 156: '(n,4np)', 157: '(n,3nd)',
    158: '(n,nda)', 159: '(n,2npa)', 160: '(n,7n)', 161: '(n,8n)',
    162: '(n,5np)', 163: '(n,6np)', 164: '(n,7np)', 165: '(n,4na)',
    166: '(n,5na)', 167: '(n,6na)', 168: '(n,7na)', 169: '(n,4nd)',
    170: '(n,5nd)', 171: '(n,6nd)', 172: '(n,3nt)', 173: '(n,4nt)',
    174: '(n,5nt)', 175: '(n,6nt)', 176: '(n,2n3He)',
    177: '(n,3n3He)', 178: '(n,4n3He)', 179: '(n,3n2p)',
    180: '(n,3n2a)', 181: '(n,3npa)', 182: '(n,dt)',
    183: '(n,npd)', 184: '(n,npt)', 185: '(n,ndt)',
    186: '(n,np3He)', 187: '(n,nd3He)', 188: '(n,nt3He)',
    189: '(n,nta)', 190: '(n,2n2p)', 191: '(n,p3He)',
    192: '(n,d3He)', 193: '(n,3Hea)', 194: '(n,4n2p)',
    195: '(n,4n2a)', 196: '(n,4npa)', 197: '(n,3p)',
    198: '(n,n3p)', 199: '(n,3n2pa)', 200: '(n,5n2p)', 203: '(n,Xp)',
    204: '(n,Xd)', 205: '(n,Xt)', 206: '(n,X3He)', 207: '(n,Xa)',
    301: 'heating', 444: 'damage-energy',
    649: '(n,pc)', 699: '(n,dc)', 749: '(n,tc)', 799: '(n,3Hec)',
    849: '(n,ac)', 891: '(n,2nc)'
}
REACTION_NAME.update({i: f'(n,n{i - 50})' for i in range(51, 91)})
REACTION_NAME.update({i: f'(n,p{i - 600})' for i in range(600, 649)})
REACTION_NAME.update({i: f'(n,d{i - 650})' for i in range(650, 699)})
REACTION_NAME.update({i: f'(n,t{i - 700})' for i in range(700, 749)})
REACTION_NAME.update({i: f'(n,3He{i - 750})' for i in range(750, 799)})
REACTION_NAME.update({i: f'(n,a{i - 800})' for i in range(800, 849)})
REACTION_NAME.update({i: f'(n,2n{i - 875})' for i in range(875, 891)})

REACTION_MT = {name: mt for mt, name in REACTION_NAME.items()}
REACTION_MT['total'] = 1
REACTION_MT['elastic'] = 2
REACTION_MT['fission'] = 18
REACTION_MT['absorption'] = 27
REACTION_MT['capture'] = 102

FISSION_MTS = (18, 19, 20, 21, 38)


def _get_products(material: Material, MT: int) -> List[Product]:
    """Generate products from MF=6 in an ENDF evaluation

    Parameters
    ----------
    material
        ENDF material to read from
    MT
        The MT value of the reaction to get products for

    Returns
    -------
    Products of the reaction

    """
    product_data = material[6, MT]
    products = []
    for data in product_data['products']:
        # Determine name of particle
        ZA = data['ZAP']
        if ZA == 0:
            name = 'photon'
        elif ZA == 1:
            name = 'neutron'
        elif ZA == 1000:
            name = 'electron'
        else:
            Z, A = divmod(ZA, 1000)
            name = gnds_name(Z, A)

        y_i = data['y_i']

        # TODO: Read distributions
        LAW = data['LAW']

        products.append(Product(name, y_i))

    return products


def _get_fission_products_endf(material: Material, MT: int) -> Tuple[List[Product], List[Product]]:
    """Generate fission products from an ENDF evaluation

    Parameters
    ----------
    ev : openmc.data.endf.Evaluation

    Returns
    -------
    products : list of openmc.data.Product
        Prompt and delayed fission neutrons
    derived_products : list of openmc.data.Product
        "Total" fission neutron

    """
    products = []
    derived_products = []

    if (1, 456) in material:
        # Prompt nu values
        data = material[1, 456]
        LNU = data['LNU']
        if LNU == 1:
            # Polynomial representation
            yield_ = Polynomial(data['C'])
        elif LNU == 2:
            # Tabulated representation
            yield_ = data['nu']

        prompt_neutron = Product('neutron', yield_=yield_)
        products.append(prompt_neutron)

    if (1, 452) in material:
        # Total nu values
        data = material[1, 452]
        if data['LNU'] == 1:
            # Polynomial representation
            yield_ = Polynomial(data['C'])
        elif data['LNU'] == 2:
            # Tabulated representation
            yield_ = data['nu']

        total_neutron = Product('neutron', yield_=yield_)
        total_neutron.emission_mode = 'total'

        if (1, 456) in material:
            derived_products.append(total_neutron)
        else:
            products.append(total_neutron)

    if (1, 455) in material:
        data = material[1, 455]
        if data['LDG'] == 0:
            # Delayed-group constants energy independent
            decay_constants = data['lambda']
            for constant in data['lambda']:
                delayed_neutron = Product('neutron')
                delayed_neutron.emission_mode = 'delayed'
                delayed_neutron.decay_rate = constant
                products.append(delayed_neutron)
        elif data['LDG'] == 1:
            # Delayed-group constants energy dependent
            raise NotImplementedError('Delayed neutron with energy-dependent '
                                      'group constants.')

        # In MF=1, MT=455, the delayed-group abundances are actually not
        # specified if the group constants are energy-independent. In this case,
        # the abundances must be inferred from MF=5, MT=455 where multiple
        # energy distributions are given.
        if data['LNU'] == 1:
            # Nu represented as polynomial
            for neutron in products[-6:]:
                neutron.yield_ = Polynomial(data['C'])
        elif data['LNU'] == 2:
            # Nu represented by tabulation
            for neutron in products[-6:]:
                neutron.yield_ = deepcopy(data['nu'])

        if (5, 455) in material:
            mf5_data = material[5, 455]
            NK = mf5_data['NK']
            if NK > 1 and len(decay_constants) == 1:
                # If only one precursor group is listed in MF=1, MT=455, use the
                # energy spectra from MF=5 to split them into different groups
                for _ in range(NK - 1):
                    products.append(deepcopy(products[1]))
            elif NK != len(decay_constants):
                raise ValueError(
                    'Number of delayed neutron fission spectra ({}) does not '
                    'match number of delayed neutron precursors ({}).'.format(
                        NK, len(decay_constants)))
            for i, subsection in enumerate(mf5_data['subsections']):
                dist = UncorrelatedAngleEnergy()
                dist.energy = EnergyDistribution.from_dict(subsection)

                delayed_neutron = products[1 + i]
                yield_ = delayed_neutron.yield_

                # Here we handle the fact that the delayed neutron yield is the
                # product of the total delayed neutron yield and the
                # "applicability" of the energy distribution law in file 5.
                applicability = subsection['p']
                if isinstance(yield_, Tabulated1D):
                    if np.all(applicability.y == applicability.y[0]):
                        yield_.y *= applicability.y[0]
                    else:
                        # Get union energy grid and ensure energies are within
                        # interpolable range of both functions
                        max_energy = min(yield_.x[-1], applicability.x[-1])
                        energy = np.union1d(yield_.x, applicability.x)
                        energy = energy[energy <= max_energy]

                        # Calculate group yield
                        group_yield = yield_(energy) * applicability(energy)
                        delayed_neutron.yield_ = Tabulated1D(energy, group_yield)
                elif isinstance(yield_, Polynomial):
                    if len(yield_) == 1:
                        delayed_neutron.yield_ = deepcopy(applicability)
                        delayed_neutron.yield_.y *= yield_.coef[0]
                    else:
                        if np.all(applicability.y == applicability.y[0]):
                            yield_.coef[0] *= applicability.y[0]
                        else:
                            raise NotImplementedError(
                                'Total delayed neutron yield and delayed group '
                                'probability are both energy-dependent.')

                delayed_neutron.distribution.append(dist)

    return products, derived_products


def _get_activation_products(material: Material, MT: int, xs: Tabulated1D) -> List[Product]:
    """Generate activation products from an ENDF evaluation

    Parameters
    ----------
    material
        ENDF material to read from
    MT
        The MT value of the reaction to get products for
    xs
        Cross section for the reaction

    Returns
    -------
    Activation products

    """
    # Determine if file 9/10 are present
    data = material[8, MT]
    present = {9: False, 10: False}
    for subsection in data['subsections']:
        if subsection['LMF'] == 9:
            present[9] = True
        elif subsection['LMF'] == 10:
            present[10] = True

    products = []

    for MF in (9, 10):
        if not present[MF]:
            continue

        data = material[MF, MT]
        for level in data['levels']:
            # Determine what the product is
            Z, A = divmod(level['IZAP'], 1000)
            excited_state = level['LFS']

            # Get GNDS name for product
            symbol = ATOMIC_SYMBOL[Z]
            if excited_state > 0:
                name = f'{symbol}{A}_e{excited_state}'
            else:
                name = f'{symbol}{A}'

            if MF == 9:
                yield_ = level['Y']
            else:
                # Re-interpolate production cross section and neutron cross
                # section to union energy grid
                production_xs = level['sigma']
                energy = np.union1d(production_xs.x, xs.x)
                prod_xs = production_xs(energy)
                neutron_xs = xs(energy)
                idx = np.where(neutron_xs > 0)

                # Calculate yield as ratio
                yield_ = np.zeros_like(energy)
                yield_[idx] = prod_xs[idx] / neutron_xs[idx]
                yield_ = Tabulated1D(energy, yield_)

            p = Product(name, yield_)
            products.append(p)

    return products


class Reaction:
    """A nuclear reaction

    This class represents a single reaction channel for a nuclide with
    an associated cross section and, if present, a secondary angle and energy
    distribution.

    Parameters
    ----------
    mt : int
        The ENDF MT number for this reaction.
    xs : dict
        Microscopic cross section for this reaction as a function of incident
        energy; these cross sections are provided in a dictionary where the key
        is the temperature of the cross section set.
    products : list
        Reaction products
    q_reaction : float
        The reaction Q-value in [eV].
    q_massdiff : float
        The mass-difference Q value in [eV].
    redundant : bool
        Indicates whether or not this is a redundant reaction

    Attributes
    ----------
    MT : int
        The ENDF MT number for this reaction.
    products : list
        Reaction products
    q_reaction : float
        The reaction Q-value in [eV].
    q_massdiff : float
        The mass-difference Q value in [eV].
    redundant : bool
        Indicates whether or not this is a redundant reaction
    xs : dict
        Microscopic cross section for this reaction as a function of incident
        energy; these cross sections are provided in a dictionary where the key
        is the temperature of the cross section set.

    """
    def __init__(self, MT: int, xs: dict = None, products: List[Product] = None,
                 q_reaction: float = 0.0, q_massdiff: float = 0.0,
                 redundant: bool = False):
        self.MT = MT
        self.xs = xs
        self.products = products
        self.q_reaction = q_reaction
        self.q_massdiff = q_massdiff
        self.redundant = redundant

    @classmethod
    def from_endf(cls, MT: int, material: Material) -> Reaction:
        """Generate reaction from ENDF file

        Parameters
        ----------
        MT
            MT value of the reaction
        material
            ENDF

        """
        # Get reaction cross section and Q values from MF=3
        rx = material[3, MT]
        q_massdiff = rx['QM']
        q_reaction = rx['QI']
        # TODO: Do something with breakup reaction flag
        xs = {'0K': rx['sigma']}

        # Get fission product yields (nu) as well as delayed neutron energy
        # distributions
        products = []
        if MT in FISSION_MTS:
            products, derived_products = _get_fission_products_endf(material, MT)
            # TODO: Store derived products somewhere

        if (6, MT) in material:
            # Product angle-energy distribution
            for product in _get_products(material, MT):
                # If fission neutrons were already added from MF=1 data, copy
                # the distribution to the existing products. Otherwise, add the
                # product to the reaction.
                if MT in FISSION_MTS and product.name == 'neutron':
                    products[0].applicability = product.applicability
                    products[0].distribution = product.distribution
                else:
                    products.append(product)

        elif (4, MT) in material or (5, MT) in material:
            # Uncorrelated angle-energy distribution
            neutron = Product('neutron')

            # Note that the energy distribution for MT=455 is read in
            # _get_fission_products_endf rather than here
            if (5, MT) in material:
                data = material[5, MT]
                for subsection in data['subsections']:
                    dist = UncorrelatedAngleEnergy()
                    dist.energy = EnergyDistribution.from_dict(subsection)

                    neutron.applicability.append(subsection['p'])
                    neutron.distribution.append(dist)
            elif MT == 2:
                # Elastic scattering -- no energy distribution is given since it
                # can be calulcated analytically
                dist = UncorrelatedAngleEnergy()
                neutron.distribution.append(dist)
            elif MT >= 51 and MT < 91:
                # Level inelastic scattering -- no energy distribution is given
                # since it can be calculated analytically. Here we determine the
                # necessary parameters to create a LevelInelastic object
                dist = UncorrelatedAngleEnergy()

                A = material[1, 451]['AWR']
                threshold = (A + 1.)/A*abs(q_reaction)
                mass_ratio = (A/(A + 1.))**2
                dist.energy = LevelInelastic(threshold, mass_ratio)

                neutron.distribution.append(dist)

            if (4, MT) in material:
                data = material[4, MT]
                for dist in neutron.distribution:
                    dist.angle = AngleDistribution.from_dict(data)

            if MT in FISSION_MTS and (5, MT) in material:
                # For fission reactions,
                products[0].applicability = neutron.applicability
                products[0].distribution = neutron.distribution
            else:
                products.append(neutron)

        if (8, MT) in material:
            for act_product in _get_activation_products(material, MT, rx['sigma']):
                # Check if product already exists from MF=6 and if it does, just
                # overwrite the existing yield.
                for product in products:
                    if act_product.name == product.name:
                        product.yield_ = act_product.yield_
                        break
                else:
                    products.append(act_product)

        return cls(MT, xs, products, q_reaction, q_massdiff)

    @classmethod
    def from_ace(cls, table: ace.Table, i_reaction: int):
        """Generate incident neutron continuous-energy data from an ACE table

        Parameters
        ----------
        table
            ACE table to read from
        i_reaction
            Index of the reaction in the ACE table

        Returns
        -------
        Reaction data

        """
        # Get nuclide energy grid
        n_grid = table.nxs[3]
        grid = table.xss[table.jxs[1]:table.jxs[1] + n_grid]*EV_PER_MEV

        # Convert temperature to a string for indexing data
        strT = temperature_str(table.temperature)

        if i_reaction > 0:
            # Get MT value
            MT = int(table.xss[table.jxs[3] + i_reaction - 1])

            # Get Q-value of reaction
            q_reaction = table.xss[table.jxs[4] + i_reaction - 1]*EV_PER_MEV

            # ==================================================================
            # CROSS SECTION

            # Get locator for cross-section data
            loc = int(table.xss[table.jxs[6] + i_reaction - 1])

            # Determine starting index on energy grid
            threshold_idx = int(table.xss[table.jxs[7] + loc - 1]) - 1

            # Determine number of energies in reaction
            n_energy = int(table.xss[table.jxs[7] + loc])
            energy = grid[threshold_idx:threshold_idx + n_energy]

            # Read reaction cross section
            xs = table.xss[table.jxs[7] + loc + 1:table.jxs[7] + loc + 1 + n_energy]

            # For damage energy production, convert to eV
            if MT == 444:
                xs *= EV_PER_MEV

            # Warn about negative cross sections
            if np.any(xs < 0.0):
                warn(f"Negative cross sections found for {MT=} in {table.name}.")

            tabulated_xs = {strT: Tabulated1D(energy, xs)}
            rx = Reaction(MT, tabulated_xs, q_reaction=q_reaction)

            # ==================================================================
            # YIELD AND ANGLE-ENERGY DISTRIBUTION

            # TODO: Read yield and angle-energy distribution

        else:
            # Elastic scattering
            mt = 2

            # Get elastic cross section values
            elastic_xs = table.xss[table.jxs[1] + 3*n_grid:table.jxs[1] + 4*n_grid]

            # Warn about negative elastic scattering cross section
            if np.any(elastic_xs < 0.0):
                warn(f"Negative elastic scattering cross section found for {table.name}.")

            xs = {strT: Tabulated1D(grid, elastic_xs)}

            # No energy distribution for elastic scattering
            # TODO: Create product

            rx = Reaction(2, xs)

        # ======================================================================
        # ANGLE DISTRIBUTION (FOR UNCORRELATED)

        # TODO: Read angular distribution

        # ======================================================================
        # PHOTON PRODUCTION

        # TODO: Read photon production

        return rx


    def __repr__(self):
        name = REACTION_NAME.get(self.MT)
        if name is not None:
            return f"<Reaction: MT={self.MT} {name}>"
        else:
            return f"<Reaction: MT={self.MT}>"
