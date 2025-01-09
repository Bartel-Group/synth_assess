import os
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
# from solidstatesynth.core.utils import in_icsd
from pydmclab.core.query import MPRester
from pydmclab.data.thermochem import gas_thermo_data
from emmet.core.thermo import ThermoType
DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"


def get_MP_data(tm_data, remake = False, thermo_types = [ThermoType.GGA_GGA_U]):
    """
    Returns a list of dictionar5ies for materials project data with properties
    of interest (the json file is of the form {'data': [entries]} an additional entry key of 'data')
    for no thermo type filter, use thermo_types = []
    """
    if os.path.exists(os.path.join(DATA_DIR, '241119_mpdata.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '241119_mpdata.json'))
    mpr = MPRester()
    entries = []
    docs = mpr.materials.thermo.search(thermo_types=thermo_types,
                                    fields=["formula_pretty",
                                                "volume",
                                                "formation_energy_per_atom",
                                                "energy_per_atom",
                                                "energy_above_hull",
                                                "theoretical",
                                                "nsites",
                                                "material_id",
                                                "task_types"],
                                                )
    print(len(docs))
    for entry in docs:
        entry_new = {}
        entry_new['material_id'] = entry.material_id
        entry_new['formula']= entry.formula_pretty
        entry_new['energy_per_atom'] = entry.energy_per_atom
        entry_new['volume'] = entry.volume
        entry_new['formation_energy_per_atom'] = entry.formation_energy_per_atom
        entry_new['energy_above_hull'] = entry.energy_above_hull
        # entry_new['theoretical'] = entry.theoretical
        entry_new['nsites'] = entry.nsites
        entry_new['tm_precursor'] = is_textmined_precursor(entry.formula_pretty, tm_data)
        if entry_new['formation_energy_per_atom'] is not None:
            entries.append(entry_new)
    # ids = [entry['material_id'] for entry in entries]
    # docs2 = mpr.materials.summary.search(fields=["theoretical",'material_id'])
    # d2 = {entry.material_id: entry.theoretical for entry in docs2}
    # print(len(docs2))
    d2 = get_is_theoretical()
    for entry in entries:
        if entry['tm_precursor']:
            entry['theoretical'] = False
        elif entry['material_id'] in d2:
            entry['theoretical'] = d2[entry['material_id']]
        else:
            entry['theoretical'] = None
    data = {'data': entries}
    fjson = os.path.join(DATA_DIR, '241119_mpdata.json')
    mp = write_json(data, fjson)
    return mp

def get_is_theoretical(remake = False):
    mpr = MPRester()
    if os.path.exists(os.path.join(DATA_DIR, '241119_is_theoretical.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '241119_is_theoretical.json'))
    docs2 = mpr.materials.summary.search(fields=["theoretical",'material_id'])
    d2 = {entry.material_id: entry.theoretical for entry in docs2}
    fjson = os.path.join(DATA_DIR, '241119_is_theoretical.json')
    mpt = write_json(d2, fjson)
    return mpt

    # add material id as key to regrab from MP
def is_textmined_precursor(formula, tm_data):
    if CompTools(formula).clean in tm_data:
        return True
    return False

def get_mp_formulas(MP):
    """
    helper function to get the list of formulas in the materials project data 
    (for any version of materials project data)
    """
    # print(MP.keys())
    mp_formulas = [entry['formula'] for entry in MP]
    return list(set(mp_formulas))

def gd_state_formula(MP, formula, tm_data):
    """
    helper function to get the ground state formula for a given formula.
    Returns the dictionary entry for the ground state of a given formula
    
    """
    mp_data = [entry for entry in MP if entry['formula'] == formula]
    gd_state_entry = mp_data[0]
    for entry in mp_data:
        if entry['formation_energy_per_atom'] < gd_state_entry['formation_energy_per_atom']:
            gd_state_entry = entry
    return gd_state_entry

def get_gd_state_MP(data,with_theoretical, tm_data,remake = False):
    """
    Returns:
        [{thermo info from MP} for ground state polymorphs in MP]
    """
    if with_theoretical == False:
        filename = '241119_mp_experimental_gd.json'
    else:
        filename = '241119_mp_gd.json'
    if os.path.exists(os.path.join(DATA_DIR, filename)) and remake == False:
        return read_json(os.path.join(DATA_DIR, filename))
    mp_formulas = get_mp_formulas(data)
    gd_state_MP = []
    for formula in mp_formulas:
        formula_gd = gd_state_formula(data, formula, tm_data)
        if formula_gd:
            gd_state_MP.append(formula_gd)
    gd_MP = {'data': gd_state_MP}
    fjson = os.path.join(DATA_DIR, filename)
    gd_mp = write_json(gd_MP, fjson)
    return gd_mp

def get_mp_experimental(MP, remake = False):
    """
    Returns a list of MP entries where the data has key 'theoretical' set to False
    (the json file is of the form {'data': [entries]} with an additional entry key of 'data')
    """
    if os.path.exists(os.path.join(DATA_DIR, '241119_mp_experimental.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '241119_mp_experimental.json'))
    entries = []
    for entry in MP:
        if not entry['theoretical']:
            entries.append(entry)
    exp_MP = {'data': entries}
    fjson = os.path.join(DATA_DIR, '241119_mp_experimental.json')
    exp_MP = write_json(exp_MP, fjson)
    return exp_MP

# def get_useful_mp_data(MP,stability_cutoff = 0.05, n_els_max = 4):
#     """
#     Args: materials project data (either as is, including only ground data, and/or including only 
#     experimental data)
#     Returns: a dictionary of materials project by data filtered by energy above the hull and
#     maximum number of elements (for efficiency)
#     """
#     MP_stability= {formula: MP[formula] for formula in MP if MP[formula]['energy_above_hull'] < stability_cutoff}
#     MP_els_max = {formula: MP_stability[formula] for formula in MP_stability if len(CompTools(formula).els) <= n_els_max}
#     return MP_els_max

def get_gases_data(remake=False):
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
        # new_gas[j]: {k: (g[j][str(k)]) / (96.485 * CompTools(j).n_atoms) for k in temperatures}
        for j in gasses
    }
    # reformatting the gasses dictionary to match the format of the materials project data
    return g_new


def main():
    tm_precursors = read_json(DATA_DIR + '/textmined_precursors.json')['precursors']
    MPt = get_is_theoretical()
    print('mp theoretical done')
    MP =  get_MP_data(tm_precursors, remake=True)['data']
    print('mp done')
    MP_exp = get_mp_experimental(MP,remake=True)['data']
    print('mp exp done')
    gd_MP = get_gd_state_MP(MP,with_theoretical = True, tm_data= tm_precursors, remake=True)['data']
    print('gd_mp done')
    gd_MP_exp = get_gd_state_MP(MP_exp,with_theoretical = False, tm_data = tm_precursors, remake=True)['data']
    print('gd_mp_exp done')
    # use only experimental data to identify the ground state so that we idenitfy the
    # experimental ground state and don't miss experimentally-observed formulas
    return MPt, MP, MP_exp, gd_MP, gd_MP_exp


if __name__ == '__main__':
    MPt, MP, MP_exp,gd_MP, gd_MP_exp = main()
