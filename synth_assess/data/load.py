import os
from pydmclab.utils.handy import read_json
from pydmclab.data.thermochem import gas_thermo_data



this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "data")
results_DATA_PATH = os.path.join(this_dir, "results_data")

def mp_data():
    return read_json(os.path.join(DATA_PATH,'mp_solids_data.json'))

def mp_data_with_theoretical():
    return read_json(os.path.join(DATA_PATH, 'mp_data_with_theoretical.json'))

def tm_precursors():
    f = read_json(os.path.join(DATA_PATH,'tm_precursors.json'))
    return f['data']

def tm_rxns():
    f = read_json(os.path.join(DATA_PATH,'tm_rxns.json'))
    return f['data']

### results (use for plotting)
def tm_entries():
    f = read_json(os.path.join(results_DATA_PATH,'tm_entries.json'))
    return f['data']

def gen_data():
    f = read_json(os.path.join(results_DATA_PATH,'gen_mat_pred_gamma.json'))
    return f

def tm_rxns_with_gamma():
    f = read_json(os.path.join(results_DATA_PATH,'tm_rxns_with_gamma.json'))
    return f

def gas_data():
    """
    Returns: dictionary of the form {temperature: {gas: dG in eV/atom}}}
    from pre-existing gas thermo data from pydmclab. Note that the gas thermo
    data provides the most accurate energies for gasses. Some gas data can be
    found in materials project as well as in gas thermo data, so any formula
    must be checked to see if it is in this dictionary before checking the
    materials project data (see get_dGf_from_source below).
    NOTE that for functionality with the ReactionNetwork code, CO2 and H2O must be
    represented as such rather than as C1O2 and H2O1 (clean formulas). Additionally,
    temperature keys are floats rather than strings
    """
    g = gas_thermo_data()
    gasses = ["H2O1", "C1O2"]
    new_gas = {'H2O1': 'H2O', 'C1O2': 'CO2'}
    temperatures = [int(key) for key in g["C1O2"]]
    g_new = {
        new_gas[j]: {k: (g[j][str(k)]) / (96.485) for k in temperatures}
        for j in gasses
    }
    # reformatting the gasses dictionary to match the format of the materials project data
    return g_new