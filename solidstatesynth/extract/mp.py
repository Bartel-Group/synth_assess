import os
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
# from solidstatesynth.core.utils import in_icsd
from pydmclab.core.query import MPRester
from pydmclab.data.thermochem import gas_thermo_data
DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"


def get_MP_data(remake = False):
    """
    Returns a dicionary of materials project data with formula as key and properties
    of interest as values (with an additional entry key of 'data')
    """
    if os.path.exists(os.path.join(DATA_DIR, '240925_mpdata.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '240925_mpdata.json'))
    mpr = MPRester()
    entries = []
    docs = mpr.materials.summary.search(fields=["formula_pretty",
                                                    "volume",
                                                    "formation_energy_per_atom",
                                                    "energy_per_atom",
                                                    "energy_above_hull",
                                                    "theoretical",
                                                    "nsites"])
    print(len(docs))
    for entry in docs:
        entry_new = {}
        entry_new['formula']= entry.formula_pretty
        entry_new['energy_per_atom'] = entry.energy_per_atom
        entry_new['volume'] = entry.volume
        entry_new['formation_energy_per_atom'] = entry.formation_energy_per_atom
        entry_new['energy_above_hull'] = entry.energy_above_hull
        entry_new['theoretical'] = entry.theoretical
        entry_new['nsites'] = entry.nsites
        if entry_new['formation_energy_per_atom'] is not None:
            entries.append(entry_new)
    data = {'data': entries}
    fjson = os.path.join(DATA_DIR, '240925_mpdata.json')
    return write_json(data, fjson)

def get_mp_formulas(MP):
    """
    helper function to get the list of formulas in the materials project data 
    (for any version of materials project data)
    """
    mp_formulas = list(set([entry['formula'] for entry in MP['data']]))
    return mp_formulas

def gd_state_formula(MP, formula):
    """
    helper function to get the ground state formula for a given formula
    
    """
    mp_data = [entry for entry in MP['data'] if entry['formula'] == formula]
    gd_state_entry = mp_data[0]
    for entry in mp_data:
        if entry['formation_energy_per_atom'] < gd_state_entry['formation_energy_per_atom']:
            gd_state_entry = entry
    return gd_state_entry

def get_gd_state_MP(MP,remake = False):
    """
    Returns:
        {formula (str) :
            {thermo info from MP}} for ground state polymorphs in MP
    """
    if os.path.exists(os.path.join(DATA_DIR, '240925_mp_ground_data.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '240925_mp_ground_data.json'))
    mp_formulas = get_mp_formulas(MP)
    gd_state_MP = {}
    for formula in mp_formulas:
        formula_gd = gd_state_formula(MP, formula)
        if formula_gd:
            gd_state_MP[CompTools(formula).clean]=formula_gd
    gd_MP = {'data': gd_state_MP}
    fjson = os.path.join(DATA_DIR, '240925_mp_ground_data.json')
    gd_mp = write_json(gd_MP, fjson)
    return gd_mp['data']

def get_mp_experimental(gd_MP, remake = False):
    if os.path.exists(os.path.join(DATA_DIR, '240926_mp_experimental.json')) and remake == False:
        return read_json(os.path.join(DATA_DIR, '240926_mp_experimental.json'))
    entries = {}
    for formula in gd_MP:
        if not gd_MP[formula]['theoretical']:
            entries[formula] = gd_MP[formula]
    exp_MP = {'data': entries}
    fjson = os.path.join(DATA_DIR, '240926_mp_experimental.json')
    exp_MP = write_json(exp_MP, fjson)
    return exp_MP['data']

def get_useful_mp_data(MP,stability_cutoff = 0.05, n_els_max = 4):
    """
    Args: materials project data (either as is, including only ground data, and/or including only 
    experimental data)
    Returns: a dictionary of materials project by data filtered by energy above the hull and
    maximum number of elements (for efficiency)
    """
    MP_stability= {formula: MP[formula] for formula in MP if MP[formula]['energy_above_hull'] < stability_cutoff}
    MP_els_max = {formula: MP_stability[formula] for formula in MP_stability if len(CompTools(formula).els) <= n_els_max}
    return MP_els_max

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
        new_gas[j]: {k: (g[j][str(k)]) / (96.485 * CompTools(j).n_atoms) for k in temperatures}
        for j in gasses
    }
    # reformatting the gasses dictionary to match the format of the materials project data
    return g_new


def main():
   MP =  get_MP_data(remake=False)['data']
   gd_MP = get_gd_state_MP(MP,remake=False)['data']
   exp_MP = get_mp_experimental(gd_MP,remake=False)['data']
   return MP, gd_MP, exp_MP

if __name__ == '__main__':
    MP, gd_MP, exp_MP = main()