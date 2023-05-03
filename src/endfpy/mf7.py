from typing import TextIO

from .records import get_head_record, get_tab1_record, get_list_record, \
    get_tab2_record


def parse_mf7_mt2(file_obj: TextIO) -> dict:
    """Generate thermal scattering data from an ENDF file

    Parameters
    ----------
    ev_or_filename : openmc.data.endf.Evaluation or str
        ENDF evaluation to read from. If given as a string, it is assumed to
        be the filename for the ENDF file.

    Returns
    -------
    openmc.data.ThermalScattering
        Thermal scattering data

    """

    # Read coherent/incoherent elastic data

    # Define helper functions to avoid duplication
    def get_coherent_elastic(file_obj):
        # Get structure factor at first temperature
        params, S = get_tab1_record(file_obj)
        T, _, LT, *_ = params
        temp_data = [{'T': T, 'LT': LT, 'S': S}]

        # Get structure factor for subsequent temperatures
        for _ in range(LT):
            params, S = get_list_record(file_obj)
            T, _, LI, *_ = params
            temp_data.append({'T': T, 'LI': LI, 'S': S})

        return temp_data

    def get_incoherent_elastic(file_obj):
        params, W = get_tab1_record(file_obj)
        return {'SB': params[0], 'W': W}

    ZA, AWR, LHTR, *_ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LHTR': LHTR}

    if LHTR == 1:
        # coherent elastic
        data['coherent'] = get_coherent_elastic(file_obj)
    elif LHTR == 2:
        # incoherent elastic
        data['incoherent'] = get_incoherent_elastic(file_obj)
    elif LHTR == 3:
        # mixed coherent / incoherent elastic
        data['coherent'] = get_coherent_elastic(file_obj)
        data['incoherent'] = get_incoherent_elastic(file_obj)

    return data


def parse_mf7_mt4(file_obj: TextIO) -> dict:
    # Read incoherent inelastic data
    ZA, AWR, _, LAT, LASYM, _ = get_head_record(file_obj)
    data = {'ZA': ZA, 'AWR': AWR, 'LAT': LAT, 'LASYM': LASYM}

    # Get information about principal atom
    params, B = get_list_record(file_obj)
    _, _, LLN, _, NI, NS = params
    data['LLN'] = LLN
    data['NI'] = NI
    data['NS'] = NS
    data['B'] = B

    # Get S(alpha,beta,T)
    kTs = []
    if data['B'][0] > 0.0:
        params, data['beta_int'] = get_tab2_record(file_obj)
        data['NB'] = NB = params[5]
        data['beta_data'] = []
        for i in range(NB):
            params, S = get_tab1_record(file_obj)
            T, beta, LT, *_ = params
            temp_data = [{'T': T, 'beta': beta, 'LT': LT, 'S': S}]
            for _ in range(LT):
                params, S = get_list_record(file_obj)
                T, beta, LI, *_ = params
                temp_data.append({'T': T, 'beta': beta, 'LT': LT, 'S': S})
            data['beta_data'].append(temp_data)

    # Get effective temperature for each atom
    _, Teff = get_tab1_record(file_obj)
    data['Teff'] = [Teff]
    for i in range(NS):
        if data['B'][6*(i + 1)] == 0.0:
            _, Teff = get_tab1_record(file_obj)
            data['Teff'].append(Teff)

    return data
