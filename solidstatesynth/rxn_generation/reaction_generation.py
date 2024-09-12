from pydmclab.core.energies import ReactionEnergy, ChemPots
from pydmclab.data.thermochem import gas_thermo_data
from pydmclab.core.comp import CompTools
from pymatgen.ext.matproj import MPRester
from pymatgen.analysis import interface_reactions
from pymatgen.analysis.phase_diagram import PhaseDiagram
# from rxn_network.costs import calculators
from solidstatesynth.utils import *
from solidstatesynth.competitions import *

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"

def get_textmined_data():
    fjson = os.path.join(DATA_DIR, "240701_corrections.json")
    new_db = read_json(fjson)['data']
    new_with_data = [entry for entry in new_db if 'common' in entry]
    return new_with_data

def get_els(formula_list):
    all_els = {}
    for formula in formula_list:
        els = CompTools(formula).els
        for el in els:
            if el not in all_els:
                all_els[el] = 1
            else :
                all_els[el] += 1
    return all_els    

def get_icsd_formulas(remake = False):
    fjson = os.path.join(DATA_DIR, "icsd_mp_entries.json")
    icsd_formula_list = read_json(fjson)['data']
    fjson2 = os.path.join(DATA_DIR, "icsd_formulas.json")
    if os.path.exists(fjson2) and not remake:
        return read_json(fjson2)
    icsd_formulas = [CompTools(entry).clean for entry in icsd_formula_list]
    return write_json({'data':icsd_formulas},fjson2)['data']

# def get_icsd_formulas(icsd_entry_data, remake = True):
#     formulas = []
#     fjson = os.path.join(DATA_DIR, "icsd_formulas.json")
#     if os.path.exists(fjson) and not remake:
#         return read_json(fjson)
#     for entry in icsd_entry_data:
#         formula = StrucTools(entry['structure']).compact_formula
#         formulas.append(formula)
#     data = {'data': list(set(formulas))}
#     return write_json(data,fjson)['data']


def get_textmined_precursors(new_db,remake = False):
    fjson = os.path.join(DATA_DIR, "textmined_precursors.json")
    if os.path.exists(fjson) and not remake:
        return read_json(fjson)
    # cantclean = 0
    precursors = []
    for entry in new_db:
        if 'common' in entry:
            prec = entry['common']['precursors']
            for p in prec:
                try:
                    p = CompTools(p).clean 
                    if p not in precursors:
                        precursors.append(p)
                except:
                    # cantclean += 1
                    pass
    # print(cantclean)
    prec_new = []
    for p in precursors:
        prec_data = MPRester().summary.search(formula = p, fields = ['material_id','energy_above_hull'])
        if prec_data:
            energies = [d.energy_above_hull for d in prec_data]
            if min(energies) < 0.05:
                # print(p)
                prec_new.append(p)
    # precursors = [p for p in precursors if MPRester().summary.search(formula = p, fields = ['material_id','energy_above_hull'])]
    # precs = [p for p in precursors if in_icsd(p)]
    p_dict = {'precursors':prec_new}
    return write_json(p_dict,fjson)

def in_icsd(formula):
    mpr = MPRester()
    data = mpr.summary.search(formula = formula, fields=['material_id','database_IDs'])
    for entry in data:
        database_ids = entry.database_IDs.get("icsd",[])
        if database_ids:
            return True
        else:
            continue
    return False

def get_relevant_precursors(target,tm_precursors):
    elems = CompTools(target).els
    # elems.extend(['O'])
    elems = list(set(elems))
    print(elems)
    precs = []
    for p in tm_precursors:
        if all(a in elems for a in list(set(CompTools(p).els))):
            if CompTools(p).chemsys != CompTools(target).chemsys:
                precs.append(p)
    # if 'O' not in CompTools(target).els:
    #     precs.append('O2')
    return precs
#add oxygen to chemical space

def get_target_reactions(target,textmined_precursors,filtered_precs = True):
    # prec_sets = []
    # if 'P' in CompTools(target).els and 'O' in CompTools(target).els:
    #     temps = [300,900,1200]
    # else:
    temps = [300,900,1200,1500,1800]
    # temps = [300,900]
    if filtered_precs:
        relevant_precs = get_relevant_precursors(target,textmined_precursors)
    else:
        relevant_precs = None
    data_by_temp = {str(temp):[] for temp in temps} 
    for temp in temps:
        data_by_temp[str(temp)]=get_competition_data(target,precursors=relevant_precs,temperature=temp, environment = 'air',open = True)
    rxns = [d for d in data_by_temp['300']]
    useful_rxns = []
    for rxn in rxns:
        rxn_dict = get_reaction_dict_from_string(rxn['rxn'])
        if len(rxn_dict['products']) == 1:
            species = list(set(rxn_dict['reactants']+rxn_dict['products']))
            useful_rxns.append(species)
        elif len(rxn_dict['products']) == 2:
            byproduct = [p for p in rxn_dict['products'] if CompTools(p).clean != CompTools(target).clean][0]
            if CompTools(byproduct).clean in ['O2','H2O1','C1O2','O1']:
                species = list(set(rxn_dict['reactants']+rxn_dict['products']))
                useful_rxns.append(species)
    useful_data = {temp:[] for temp in data_by_temp}
    for temp in useful_data:
        for d in data_by_temp[temp]:
            rxn_dict_new = get_reaction_dict_from_string(d['rxn'])
            spec =  list(set(rxn_dict_new['reactants']+rxn_dict_new['products']))
            if spec in useful_rxns:
                useful_data[temp].append(d)
        # print(useful_data.keys())
        # useful_data[temp] = [d for d in data_by_temp[temp] if get_reaction_dict_from_string(d['rxn']) in useful_rxns]
        # print('temp',len(useful_data[temp]))
    return useful_data

def get_filtered_target_reactions(target, textmined_precursors,filtered_precs = True):
    useful_data = get_target_reactions(target,textmined_precursors,filtered_precs)
    min_gamma_data = {}
    if not useful_data['300']:
        return useful_data
    elif len(useful_data['300']) == 0:
        return useful_data
    for temp in useful_data:
        min_data = useful_data[temp][0]
        for d in useful_data[temp]:
            if d['gamma'] < min_data['gamma']:
                min_data = d
        min_gamma_data[temp] = min_data
    return min_gamma_data

def get_icsd_formulas_reduced(icsd_formulas):
    icsd_formulas_reduced = []
    for formula in icsd_formulas:
        if CompTools(formula).n_els==3:
            if all(el not in ['Tc', 'Pa', 'Xe', 'Pm', 'Kr'] for el in CompTools(formula).els):
                icsd_formulas_reduced.append(formula)
    return icsd_formulas_reduced

def get_icsd_targets_by_chemsys(icsd_formulas_reduced):
    icsd_targets = {}
    for formula in icsd_formulas_reduced:
        chemsys = CompTools(formula).chemsys
        if chemsys not in icsd_targets:
            icsd_targets[chemsys] = [formula]
        else:
            icsd_targets[chemsys].append(formula)
    return icsd_targets

def get_precursors_by_chemsys(icsd_targets_by_chemsys,tm_precursors,remake=False):
    fjson = os.path.join(DATA_DIR, "240910_chemsys_precursors.json")
    if os.path.exists(fjson) and not remake:
        return read_json(fjson)
    precs_by_chemsys = {}
    for chemsys in icsd_targets_by_chemsys:
        target = icsd_targets_by_chemsys[chemsys][0]
        relevant_precs = get_relevant_precursors(target,tm_precursors)
        precs_by_chemsys[chemsys] = relevant_precs
    data = {'data':precs_by_chemsys}
    return write_json(data,fjson)


def get_optimum_reactions(icsd_formula,textmined_precursors):
    # chemsys = CompTools(icsd_formula).chemsys
    try:
        filtered_reactions = get_filtered_target_reactions(icsd_formula,textmined_precursors)
    except:
        try:
            filtered_reactions = get_filtered_target_reactions(icsd_formula,textmined_precursors,filtered_precs = False)
            print('using all precursors',icsd_formula)        
        except:
            print('alternative error',icsd_formula)
            filtered_reactions = {}
    # except AttributeError as e:
    #     print('attribute error',icsd_formula)
    #     filtered_reactions = get_filtered_target_reactions(icsd_formula,textmined_precursors,filtered_precs = False)
    # except KeyError as e:
    #     print('key error',icsd_formula)
    #     filtered_reactions = {}
        # filtered_reactions = {}
    return {icsd_formula:filtered_reactions}

def get_optimum_reactions_for_all_icsd_targets(icsd_formulas_reduced,textmined_precursors,n_procs = 4, remake = True):
    icsd_rxns = {}
    fjson = os.path.join(DATA_DIR, "240910_icsd_reactions.json")
    if os.path.exists(fjson) and not remake:
        return read_json(fjson)
    # reaction_indices = list(range(len(reaction_data)))

    pool = mp.Pool(processes=n_procs)
    results = [
        r
        for r in pool.starmap(
            get_optimum_reactions,
            [
                (
                    formula,
                    # precs_by_chemsys,
                    textmined_precursors,
                )
                for formula in icsd_formulas_reduced
            ],
        )
    ]
    pool.close()
    for entry in results:
        formula = list(entry.keys())[0]
        icsd_rxns[formula] = entry[formula]
    icsd_rxns = {'data':icsd_rxns}
    entry_json = write_json(icsd_rxns, fjson)
    return entry_json['data']





def main():
    new_db = get_textmined_data()
    textmined_precursors = get_textmined_precursors(new_db)['precursors']
    icsd_formulas = get_icsd_formulas(remake = False)['data']
    icsd_formulas_reduced = get_icsd_formulas_reduced(icsd_formulas)
    icsd_targets_by_chemsys = get_icsd_targets_by_chemsys(icsd_formulas_reduced)
    precs_by_chemsys = get_precursors_by_chemsys(icsd_targets_by_chemsys,textmined_precursors)['data']
    icsd_with_oxygen = [f for f in icsd_formulas_reduced if 'O' in CompTools(f).els]
    # icsd_rxns = get_optimum_reactions_for_all_icsd_targets(icsd_with_oxygen[0:10],precs_by_chemsys,textmined_precursors,n_procs = 4, remake = False)

    # useful_data = get_target_reactions
    return new_db, textmined_precursors, icsd_formulas, icsd_formulas_reduced, precs_by_chemsys, icsd_with_oxygen


if __name__ == "__main__":
    (new_db,
     textmined_precursors,
     icsd_formulas, 
     icsd_formulas_reduced,
     precs_by_chemsys,
     icsd_with_oxygen,
    ) = main()     