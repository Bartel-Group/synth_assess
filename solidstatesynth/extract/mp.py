import os
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
# from solidstatesynth.core.utils import in_icsd
from pydmclab.core.query import MPRester
DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"


def get_MP_data(remake = False):
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
    mp_formulas = list(set([entry['formula'] for entry in MP['data']]))
    return mp_formulas

def gd_state_formula(MP, formula):
    """
    
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
            {thermo info from MP}}
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

def main():
   MP =  get_MP_data(remake=False)
   gd_MP = get_gd_state_MP(MP,remake=True)
   exp_MP = get_mp_experimental(gd_MP,remake=True)
   return MP, gd_MP, exp_MP

if __name__ == '__main__':
    MP, gd_MP, exp_MP = main()