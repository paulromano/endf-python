# SPDX-FileCopyrightText: 2023 Paul Romano
# SPDX-License-Identifier: MIT

from typing import TextIO

import numpy as np

from .records import get_head_record, get_cont_record, get_tab1_record, \
    get_list_record



def parse_mf2(file_obj: TextIO) -> dict:
    """Parse resonance parameters from MF=2, MT=151

    Parameters
    ----------
    file_obj
        File-like object to read from

    Returns
    -------
    dict
        Resonance parameter data

    """
    # Determine whether discrete or continuous representation
    ZA, AWR, _, _, NIS, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'NIS': NIS}

    data['isotopes'] = []
    for _ in range(NIS):
        ZAI, ABN, _, LFW, NER, _ = get_cont_record(file_obj)
        iso = {'ZAI': ZAI, 'ABN': ABN, 'LFW': LFW, 'NER': NER}

        iso['ranges'] = []
        for _ in range(NER):
            EL, EH, LRU, LRF, NRO, NAPS = get_cont_record(file_obj)
            rrange = {'EL': EL, 'EH': EH, 'LRU': LRU,
                      'LRF': LRF, 'NRO': NRO, 'NAPS': NAPS}

            if LRF == 0:
                # Read spin and scattering radius
                SPI, AP, _, _, NLS, _ = get_cont_record(file_obj)
                rrange['SPI'] = SPI
                rrange['AP'] = AP
                rrange['NLS'] = NLS
            elif LRU in (0, 1):
                # resolved resonance region
                rrange.update(_FORMALISMS[LRF].dict_from_endf(file_obj, NRO))
            elif LRF == 2:
                # unresolved resonance region
                rrange.update(Unresolved.dict_from_endf(file_obj, LFW, LRF, NRO))
            iso['ranges'].append(rrange)

        data['isotopes'].append(iso)

    return data


class MLBW:
    @staticmethod
    def dict_from_endf(file_obj: TextIO, NRO: int) -> dict:
        # Read energy-dependent scattering radius if present
        data = {}
        if NRO != 0:
            _, data['APE'] = get_tab1_record(file_obj)

        # Other scatter radius parameters
        SPI, AP, _, _, NLS, _ = get_cont_record(file_obj)
        data['SPI'] = SPI
        data['AP'] = AP
        data['NLS'] = NLS

        # Read resonance widths, J values, etc
        data['sections'] = []
        for l in range(NLS):
            (AWRI, QX, L, LRX, _, NRS), values = get_list_record(file_obj)
            section = {'AWRI': AWRI, 'QX': QX, 'L': L, 'LRX': LRX, 'NRS': NRS}
            section['ER'] = values[0::6]
            section['AJ'] = values[1::6]
            section['GT'] = values[2::6]
            section['GN'] = values[3::6]
            section['GG'] = values[4::6]
            section['GF'] = values[5::6]
            data['sections'].append(section)

        return data


class SLBW(MLBW):
    ...


class ReichMoore:
    @staticmethod
    def dict_from_endf(file_obj: TextIO, NRO: int) -> dict:
        # Read energy-dependent scattering radius if present
        data = {}
        if NRO != 0:
            _, data['APE'] = get_tab1_record(file_obj)

        # Other scatter radius parameters
        SPI, AP, LAD, _, NLS, NLSC = get_cont_record(file_obj)
        data['SPI'] = SPI
        data['AP'] = AP
        data['LAD'] = LAD
        data['NLS'] = NLS
        data['NLSC'] = NLSC

        # Read resonance widths, J values, etc
        data['sections'] = []
        for l in range(NLS):
            (AWRI, APL, L, _, _, NRS), values = get_list_record(file_obj)
            section = {'AWRI': AWRI, 'APL': APL, 'L': L, 'NRS': NRS}
            section['ER'] = values[0::6]
            section['AJ'] = values[1::6]
            section['GN'] = values[2::6]
            section['GG'] = values[3::6]
            section['GFA'] = values[4::6]
            section['GFB'] = values[5::6]
            data['sections'].append(section)

        return data


class AdlerAdler:
    @staticmethod
    def dict_from_endf(file_obj: TextIO, NRO: int) -> dict:
        raise NotImplementedError


class RMatrixLimited:
    @staticmethod
    def dict_from_endf(file_obj: TextIO, NRO: int) -> dict:
        _, _, IFG, KRM, NJS, KRL = get_cont_record(file_obj)
        data = {'IFG': IFG, 'KRM': KRM, 'NJS': NJS, 'KRL': KRL}

        # Read particle-pair data
        items, values = get_list_record(file_obj)
        data['NPP'] = NPP = items[2]
        data['particle_pairs'] = pp = {}
        pp['MA'] = values[::12]
        pp['MB'] = values[1::12]
        pp['ZA'] = values[2::12]
        pp['ZB'] = values[3::12]
        pp['IA'] = values[4::12]
        pp['IB'] = values[5::12]
        pp['Q'] = values[6::12]
        pp['PNT'] = values[7::12]
        pp['SHF'] = values[8::12]
        pp['MT'] = values[9::12]
        pp['PA'] = values[10::12]
        pp['PB'] = values[11::12]

        # loop over spin groups
        data['spin_groups'] = spin_groups = []
        for i in range(NJS):
            (AJ, PJ, KBK, KPS, _, NCH), values = get_list_record(file_obj)
            spin_group = {'AJ': AJ, 'PJ': PJ, 'KBK': KBK, 'KPS': KPS, 'NCH': NCH}
            spin_group['channels'] = channels = {}
            channels['PPI'] = values[::6]
            channels['L'] = values[1::6]
            channels['SCH'] = values[2::6]
            channels['BND'] = values[3::6]
            channels['APE'] = values[4::6]
            channels['APT'] = values[5::6]

            # Read resonance energies and widths
            (*_, NRS, _, NX), values = get_list_record(file_obj)
            spin_group['NRS'] = NRS
            spin_group['NX'] = NX
            spin_group['ER'] = values[::NCH + 1]

            # Read widths into a matrix and transpose
            GAM = []
            for j in range(NRS):
                GAM.append(values[1 + (NCH + 1)*j:(NCH + 1)*(j + 1)])
            GAM = np.array(GAM).reshape(NRS, NCH)
            spin_group['GAM'] = GAM.T

            # Optional extension (Background R-Matrix)
            if KBK > 0:
                (_, _, LCH, LBK, _, _), values = get_list_record(file_obj)
                spin_group['LCH'] = LCH
                spin_group['LBK'] = LBK
                if LBK == 1:
                    _, spin_group['RBR'] = get_tab1_record(file_obj)
                    _, spin_group['RBI'] = get_tab1_record(file_obj)
                elif LBK in (2, 3):
                    (ED, EU, *_), values = get_list_record(file_obj)
                    spin_group['ED'] = ED
                    spin_group['EU'] = EU

            # Optional extension (Tabulated phase shifts)
            if KPS > 0:
                items, values = get_list_record(file_obj)
                spin_group['LPS'] = LPS = items[4]
                if LPS == 1:
                    _, spin_group['PSR'] = get_tab1_record(file_obj)
                    _, spin_group['PSI'] = get_tab1_record(file_obj)

            spin_groups.append(spin_group)

        return data


class Unresolved:
    @staticmethod
    def dict_from_endf(file_obj: TextIO, LFW: int, LRF: int, NRO: int) -> dict:
        # Read energy-dependent scattering radius if present
        data = {}
        if NRO != 0:
            _, data['APE'] = get_tab1_record(file_obj)

        # Get SPI, AP, and LSSF
        if not (LFW == 1 and LRF == 1):
            SPI, AP, LSSF, _, NLS, _ = get_cont_record(file_obj)
            data.update({'SPI': SPI, 'AP': AP, 'LSSF': LSSF, 'NLS': NLS})

        if LFW == 0 and LRF == 1:
            # Case A -- fission widths not given, all parameters are
            # energy-independent
            data['ranges'] = []
            for _ in range(NLS):
                (AWRI, _, L, _, _, NJS), values = get_list_record(file_obj)
                rrange = {'AWRI': AWRI, 'L': L, 'NJS': NJS}
                rrange['D'] = values[::6]
                rrange['AJ'] = values[1::6]
                rrange['AMUN'] = values[2::6]
                rrange['GNO'] = values[3::6]
                rrange['GG'] = values[4::6]
                data['ranges'].append(rrange)

        elif LFW == 1 and LRF == 1:
            # Case B -- fission widths given, only fission widths are
            # energy-dependent
            (SPI, AP, LSSF, _, NE, NLS), ES = get_list_record(file_obj)
            data.update({'SPI': SPI, 'AP': AP, 'LSSF': LSSF,
                         'NE': NE, 'NLS': NLS, 'ES': ES})

            data['ranges'] = []
            for _ in range(NLS):
                AWRI, _, L, _, NJS, _ = get_cont_record(file_obj)
                rrange = {'AWRI': AWRI, 'L': L, 'NJS': NJS}
                rrange['parameters'] = []
                for j in range(NJS):
                    items, values = get_list_record(file_obj)
                    rrange['parameters'].append({
                        'MUF': items[3], 'D': values[0], 'AJ': values[1],
                        'AMUN': values[2], 'GN0': values[3], 'GG': values[4],
                        'GF': values[6:]
                    })
                data['ranges'].append(rrange)

        elif LRF == 2:
            # Case C -- all parameters are energy-dependent
            data['ranges'] = []
            for _ in range(NLS):
                AWRI, _, L, _, NJS, _ = get_cont_record(file_obj)
                rrange = {'AWRI': AWRI, 'L': L, 'NJS': NJS}
                rrange['parameters'] = []
                for _ in range(NJS):
                    (AJ, _, INT, _, _, NE), values = get_list_record(file_obj)
                    rrange['parameters'].append({
                        'AJ': AJ, 'INT': INT, 'NE': NE,
                        'AMUX': values[2],
                        'AMUN': values[3],
                        'AMUF': values[5],
                        'E': values[6::6],
                        'D': values[7::6],
                        'GX': values[8::6],
                        'GN0': values[9::6],
                        'GG': values[10::6],
                        'GF': values[11::6]
                    })
                data['ranges'].append(rrange)

        return data


_FORMALISMS = {
    1: SLBW,
    2: MLBW,
    3: ReichMoore,
    4: AdlerAdler,
    7: RMatrixLimited
}
