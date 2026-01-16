import os
from pydmclab.core.comp import CompTools
from solidstatesynth.gen.rxn_metrics import EnumerateRxns, RxnsAtNewTempEnv, AnalyzeReactionSet
DATA_DIR = "/user/specifies/data/path"
solids_data = os.path.join(DATA_DIR, "mp_solids_data.json")

def get_metrics(target, temperature, solids_data, gen_data = None, is_gen = False, restrict_to_tm = False):
    """
    Args:
        target: string
        temperature: reaction temperature of interest (K)
        solids_data: Dictionary of the solids data to use for the precursors, of the form {formula:{data for formula}}
        -- for filtering, data must include keys "theoretical" (bool), "textmined_precursor" (bool) and "energy_above_hull" (float)
        restrict_to_tm (bool): Whether to restrict the precursors to only those listed in the text-mined dataset
        gen_data (dict): dictionary of dictionaries with formulas as keys (gen_data should be only for a specific model type)
        is_gen (bool): True if metrics are being generated for material not in MP. In this case, gen_data is also required
    Returns:
        A dictionary of reactions and associatede computed metrics
    """
    if is_gen:
        gen_formula = target
    else:
        gen_formula = None
    r = EnumerateRxns(els = CompTools(target).els, open = False, temperature = 300, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula, prec_kwargs={'restrict_to_tm_precursors':restrict_to_tm}).rxns
    r = RxnsAtNewTempEnv(reaction_set = r, els = CompTools(target).els, new_temperature = temperature, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula).corrected_reactions_at_temp()
    print(list(r))
    t = AnalyzeReactionSet(target = target, reactions = r, temperature = temperature)
    return t.metrics_at_temp_env()



# def classifications(struc_file):
#   return {'PU-CGNF':, 'SynthNN':, 'PU-CGCNN':, 'SynCoTrain':}
# lines to run classification models GO HERE


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
                        CompTools(reactant).clean
                        reactant_list.append(reactant)
                        reactant_coeffs.append(coefficient)
                    except:
                        return None
        else:
            reactant_list.append(reactant)
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
                        CompTools(product).clean
                        product_list.append(product)
                        product_coeffs.append(coefficient)
                    except:
                        return None
        else:
            product_list.append(product)
            product_coeffs.append(1)
    # print(reactant_list,product_list)
    return {"reactants": reactant_list, "products": product_list, 
            "reactant_coeffs": reactant_coeffs, "product_coeffs": product_coeffs}

