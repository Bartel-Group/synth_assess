
from pydmclab.core.comp import CompTools
from pydmclab.core.query import MPRester
from pydmclab.core.energies import ReactionEnergy
from rxn_network.reactions.computed import ComputedReaction as CR
from rxn_network.entries.entry_set import GibbsEntrySet
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry,GibbsComputedStructureEntry
import os
import math
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import re


DATA_DIR_MP = "/Volumes/cems_bartel/projects/mp-stability/data"
DATA_DIR_paperdict = "/Volumes/cems_bartel/projects/negative-examples/data"
api_key = 'fDY65nUJRZpelqz70f36AHv3cxv2EpkY'

def in_icsd(formula):
    """
    Returns True if the formula is in the ICSD database
    """
    mpr = MPRester()
    print('formula',formula)
    data = mpr.materials.summary.search(formula = formula, fields=['material_id','database_IDs'])
    for entry in data:
        database_ids = entry.database_IDs.get("icsd",[])
        if database_ids:
            return True
        else:
            continue
    return False

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

def get_chemical_potential(env,temp):
    std_pres = 1E5
    boltzmann_const = 8.617E-5

    if env in ['nitrogen', 'argon', 'hydrogen', 'N2','Ar','H2', 'inert', 
                 'carbon monoxide', 'CO','carbon','platinum',]:
        pp_dict = {'O': 1E-6,
               'C1O2':1E-6,
               'H2O1':1E-6}
    elif env in ['oxygen','O2']:
        pp_dict = {'O': std_pres,
               'C1O2':1E-6,
               'H2O1':1E-6
               }
    else:
        pp_dict = {'O': 2E4,
               'C1O2':40,
               'H2O1':2.3E3
               }
    pres = pp_dict['O']
    mu = boltzmann_const * temp * np.log(pres / std_pres)
    return mu

# def get_ComputedEntry(formula):

def get_gibbs_entry_set(elems,temperature):
    mpr = MPRester(api_key)
    mp_entries = mpr.get_entries_in_chemsys(elems)
    gibbs_entry_set = GibbsEntrySet.from_computed_entries(mp_entries, temperature,include_nist_data=True)
    gibbs_entry_set = gibbs_entry_set.filter_by_stability(e_above_hull=0.05)
    return gibbs_entry_set

def get_balanced_reaction_coefficients(precursors, target):
    formulas = precursors + [target]
    dGf_dict = {str(formula): {"E": 0} for formula in formulas}
    formulas = precursors + [target] + ['C1O2','H2O1','O2']
    dGf_dict = {str(formula): {"E": 0} for formula in formulas}
    re = ReactionEnergy(
    input_energies=dGf_dict,
    reactants=precursors,
    products=[target,'C1O2','O2','H2O1'],
    norm="atom",
    energy_key="E",
    )
    try:
        re.coefs
        return re.coefs
    except:
        # print('cannot balance',precursors,target)
        return None


def main():
    return

if __name__ == "__main__":
    main()