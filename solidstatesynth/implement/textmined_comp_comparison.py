# from solidstatesynth.gen.build_rxn import BuildRxn, AnalyzeRxnString, AnalyzeEnumeratedRxns
from solidstatesynth.extract.tm import get_updated_textmined_data,get_mp_cmpds,get_tm_rxns
from pydmclab.core.comp import CompTools
import os
from pydmclab.utils.handy import read_json, write_json
# from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeChemsys
from solidstatesynth.gen.rxn_metrics import PrecursorSet, EnumerateRxns, TempEnvCorrections, RxnsAtNewTempEnv, AnalyzeReactionSet

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"
common_chemsystems = ['Bi-Fe-O', 'Li-O-Ti', 'Li-Mn-O', 'Co-Li-O', 'O-Sr-Ti', 'Ba-O-Ti', 'Al-O-Y']

def tm_rxns():
    """
    Get all thermomechanical reactions from the Materials Project.
    """
    tm = get_updated_textmined_data()
    cmpds = get_mp_cmpds()
    tm_rxns = get_tm_rxns(tm, cmpds)
    return tm_rxns

def relevant_tm_rxns():
    """
    Get all thermomechanical reactions from the Materials Project that are relevant to the negative example generation.
    """
    rxns = tm_rxns()
    relevant_tm_rxns = []
    for rxn in rxns:
        target = rxn['target']
        target_els = CompTools(target).els
        if len(target_els)==3:
            if 'O' in target_els:
                environment = rxn['atmosphere']
                if environment == 'air' or environment == None:
                    relevant_tm_rxns.append(rxn)
    return relevant_tm_rxns

def relevant_tm_new():
    rxns = tm_rxns()
    relevant_tm_rxns = []
    for rxn in rxns:
        target = rxn['target']
        target_els = CompTools(target).els
        if 'O' in target_els:
            if all([a not in target_els for a in ['C', 'H', 'N']]):
                environment = rxn['atmosphere']
                if environment == 'air' or environment == None:
                    if len(target_els)==4:
                        relevant_tm_rxns.append(rxn)
    return relevant_tm_rxns

def chemsys_reactions(els):
    reactions = EnumerateRxns(els).rxns
    return reactions

def get_reaction_dict_from_string(reaction_string):
    """
    Args: reaction string
    Returns: dictionary with keys 'reactants' and 'products' where
    the values are lists of the reactants and products in the reaction
    *** IF CLEANABLE -- otherwise returns None ***
    Uses: this is the most usable reaction format from which to calculate
    dG_rxn using the pydmclab ReactionEnergy class-- get_dGrxn_at_T takes
    a reaction dictionary as an argument
    """
    reactant_list = []
    reactant_coeffs = []
    product_list = []
    product_coeffs = []
    if '->' in reaction_string:
        reactants, products = reaction_string.split(" -> ")
        # print(reactants)
    elif '==' in reaction_string:
        reactants, products = reaction_string.split(" == ")
    if "+" in reactants:
        reactants = reactants.split(" + ")
    else:
        reactants = [reactants]
    # print('reactants', reactants)
    for reactant in reactants:
        if " " in reactant:
            if len(reactant.split(" ")) != 1:
                coefficient, reactant = reactant.split(" ")
                if float(coefficient) > 0:
                    try:
                        reactant = CompTools(reactant).clean
                        reactant_list.append(reactant)
                        reactant_coeffs.append(coefficient)
                    except:
                        return None
        else:
            # print(reactant)
            reactant_list.append(CompTools(reactant).clean)
            reactant_coeffs.append(1)
    if "+" in products:
        products = products.split(" + ")
    else:
        products = [products]
    for product in products:
        if " " in product:
            if len(product.split(" ")) != 1:
                coefficient, product = product.split(" ")
                if float(coefficient) > 0:
                    try:
                        product = CompTools(product).clean
                        product_list.append(product)
                        product_coeffs.append(coefficient)
                    except:
                        return None
        else:
            # print(product)
            product_list.append(CompTools(product).clean)
            product_coeffs.append(1)
    # print(reactant_list,product_list)
    return {"reactants": reactant_list, "products": product_list, 
            "reactant_coeffs": reactant_coeffs, "product_coeffs": product_coeffs}


def metrics_for_tm_entry(tm_entry, solids_data, temperature_new = None):
    tm_target = tm_entry['target']
    rxn_els = CompTools(tm_target).els
    tm_rxn = tm_entry['rxn']
    rxn_dict = get_reaction_dict_from_string(tm_rxn)
    tm_precursors = rxn_dict['reactants']
    tm_environment = tm_entry['atmosphere']
    if not tm_environment:
        tm_environment = 'air'
    tm_temperature = tm_entry['temperature']
    if not tm_temperature:
        tm_temperature = 1073
    tm_rxns = EnumerateRxns(rxn_els,solids_data=solids_data, temperature = tm_temperature).rxns
    if temperature_new:
        tm_env_rxns = RxnsAtNewTempEnv(tm_rxns, temperature_new, tm_environment).corrected_reactions_at_temp
    else:
        tm_rxns_env = TempEnvCorrections(tm_temperature, tm_environment).rxns_with_temp_env_correction(tm_rxns)
    metrics = AnalyzeReactionSet(tm_rxns_env, target = tm_target, temperature=tm_temperature).metrics_at_temp_env()
    return metrics
    
def matched_rxns(tm_entry, metrics):
    tm_rxn = tm_entry['rxn']
    rxn_dict_tm = get_reaction_dict_from_string(tm_rxn)
    species_tm = set(rxn_dict_tm['reactants'] + rxn_dict_tm['products'])
    print('tm',species_tm)
    for entry in metrics:
        rxn_dict = get_reaction_dict_from_string(entry['rxn'])
        species = set(rxn_dict['reactants'] + rxn_dict['products'])
        print(species)
        if species_tm == species:
            return entry    
    return None 

def optimum_rxn(metrics):
    min_gamma = metrics[0]['gamma']
    min_entry = metrics[0]
    for entry in metrics:
        if entry['gamma'] < min_gamma:
            min_gamma = entry['gamma']
            min_entry = entry
    return min_entry

def reaction_comparisons(tm_entry, solids_data, temperature_new=None):
    metrics = metrics_for_tm_entry(tm_entry, solids_data, temperature_new)
    # print(metrics)
    matched = matched_rxns(tm_entry, metrics)
    optimum = optimum_rxn(metrics)
    return matched, optimum

def tm_rxns_by_chemsys(tm):
    rxns_by_chemsys = {}
    new = True if len(CompTools(tm[0]['target']).els) > 3 else False
    for entry in tm:
        chemsys = CompTools(entry['target']).chemsys
        chemsys = chemsys.replace('-',',')
        if chemsys not in rxns_by_chemsys:
            rxns_by_chemsys[chemsys] = []
        rxns_by_chemsys[chemsys].append(entry)
    if new:
        d = os.path.join(DATA_DIR, 'tm_chemsys_new.json')
    else:
        d = os.path.join(DATA_DIR, 'tm_chemsys.json')

    # r = write_json(rxns_by_chemsys, d)
    return write_json(rxns_by_chemsys, d)

def reduced_rxns_by_chemsys(tm):
    reduced = {}
    for chemsys in tm_rxns_by_chemsys:
        reduced[chemsys] = {}
        for entry in tm_rxns_by_chemsys[chemsys]:
            reduced[chemsys].append({'target':entry['target'], 'temperature': entry['temperature']})
    return reduced

# def get_competitions(rxns):
#     competitions = []
#     for tm_rxn in rxns:
#         optimum, actual_rxn = competitions_for_tm_entry(tm_rxn)
#         if actual_rxn:
#             competitions.append({'optimum':optimum, 'actual':actual_rxn})
#     return competitions

def main():
    tm = relevant_tm_rxns()
    tm_new = relevant_tm_new()
    fjson = os.path.join(DATA_DIR, "241119_mp_gd.json")
    solids_data = read_json(fjson)['data']
    solids_data = {CompTools(entry['formula']).clean: entry for entry in solids_data}
    tm_chemsys = tm_rxns_by_chemsys(tm)
    tm_chemsys_new = tm_rxns_by_chemsys(tm_new)
    tm_chemsys_all = tm_chemsys|tm_chemsys_new
    tm_chemsys_all = {key: [a for a in tm_chemsys_all[key] if a['mp']] for key in tm_chemsys_all}
    tm_chemsys_all = {key: tm_chemsys_all[key] for key in tm_chemsys_all if tm_chemsys_all[key]}
    f = write_json(tm_chemsys_all,  os.path.join(DATA_DIR, 'tm_chemsys_all.json'))
    return tm, solids_data, tm_chemsys, tm_chemsys_new, tm_chemsys_all

if __name__ == "__main__":
    tm, solids_data, tm_chemsys, tm_chemsys_new, tm_chemsys_all = main()