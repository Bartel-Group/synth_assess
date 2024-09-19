from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
from pydmclab.data.thermochem import gas_thermo_data
from pydmclab.core.energies import ReactionEnergy
from pydmclab.core.query import MPRester
from pydmclab.hpc.helpers import get_query
from rxn_network.reactions.hull import InterfaceReactionHull as IRH
from rxn_network.reactions.computed import ComputedReaction as CR
from rxn_network.entries.entry_set import GibbsEntrySet
from pymatgen.entries.computed_entries import ComputedEntry, ComputedStructureEntry,GibbsComputedStructureEntry
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
# from pymatgen.analysis.phase_diagram import PhaseDiagram
# from pymatgen.analysis import interface_reactions
# from pymatgen.entries import computed_entries
# from rxn_network.costs import calculators
# from mp_api.client import MPRester
import os
import math
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import re


DATA_DIR_MP = "/Volumes/cems_bartel/projects/mp-stability/data"
DATA_DIR_paperdict = "/Volumes/cems_bartel/projects/negative-examples/data"
api_key = 'fDY65nUJRZpelqz70f36AHv3cxv2EpkY'

def get_textmined_data():
    """
    Function: Reads text-mined paper data json file yielding a nested dictionary
    inside the 'reactions' key.
    This is the paper data associated with Kononova et al.'s text-mined dataset
    for synthesizability modified to include Chris B.'s dG_rxn calculations.
    Use: The dictionary is used to find successfully synthesized targets and the
    associated synthesis reaction, temperature, and environment. Note that dG_rxn
    is not drawn from this database because the initial dG_rxn calculations were
    performed using the old materials project data
    """
    fjson = os.path.join(DATA_DIR_paperdict, "paperdict.json")
    pd = read_json(fjson)
    return pd['reactions']


def get_mp_data():
    """
    THIS FUNCTION CALLS NEW MP DATA
    Function: Reads materials project dataset of the form
    {temperature: {formula: {Ef: dG in eV/atom}}}
    Materials Project formula keys should be cleaned formulas
    (use pydmc CompTools(formula).clean).
    The MP energies are calculated using DFT and the 230608 dataset
    has the newly calculated MP energies.
    Use: MP energies are subsequently used to calculate dG_rxn
    """
    fjson = os.path.join(DATA_DIR_MP, "230608_mp_stability.json")
    return read_json(fjson)

def get_interpolation_data():
    """
    Reads dataset with interpolated dG values for formulas not in materials project.
    The dGf associated with a given formula is interpolated from the dGf values of
    closeby materials in Materials Project.
    The interpolation data is of the form {formula:{"formula":, "interpolation":, "dGf":{temp:}}}.
    Note that this interpolation data is modified to reflect updated dGf values in the materials project.
    """
    fjson = os.path.join(DATA_DIR_paperdict, "24_modified_interpolations.json")
    return read_json(fjson)

def get_gases_data():
    """
    Returns: dictionary of the form {temperature: {gas: {Ef: dG in eV/atom}}}
    from pre-existing gas thermo data from pydmclab. Note that the gas thermo
    data provides the most accurate energies for gasses. Some gas data can be
    found in materials project as well as in gas thermo data, so any formula
    must be checked to see if it is in this dictionary before checking the
    materials project data (see get_dGf_from_source below).
    """
    g = gas_thermo_data()
    gasses = ["H2O1", "C1O2"]
    temperatures = g["C1O2"].keys()
    g_new = {
        k: {j: {"Ef": (g[j][k]) / (96.485 * CompTools(j).n_atoms)} for j in gasses}
        for k in temperatures
    }
    # reformatting the gasses dictionary to match the format of the materials project data
    return g_new

def get_updated_textmined_data():
    """
    Returns: list of dictionaries of the form
    {"common":{"doi": ..., "temperature":..., "environment":...,
    "chemsys":..., "precursors":...}, "positive":{ "reaction": reaction,
    "target": target, "dG_rxn": dG_rxn}}
    with each item associated with an entry in "thermo" (a single reaction)
    """
    fjson = os.path.join(DATA_DIR_paperdict, "240620_nocorrections.json")
    newdb = read_json(fjson)
    return newdb['data']

def get_gas_corrected_data():

    fjson = os.path.join(DATA_DIR_paperdict, "240620_gascorrections")
    new_db = read_json(fjson)
    return new_db['data']

def get_restructured_db():

    fjson = os.path.join(DATA_DIR_paperdict, "240513_restructured_db.json")
    restructured_db = read_json(fjson)
    return restructured_db

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

def get_dGf_from_source(formula, materials_project_data, gas_data, interpolation_data):
    """
    Args: formula string (does NOT need to be pre-cleaned), databases
    Returns: dictionary of the form {temperature: dGf} or None if formula not found
    Function: This function determine which database should be used for a given formula and returns
    the dGf values for that formula at all temperatures in the database.
    NOTE THAT some formulas may be found in gasses and in materials project. In this case, gasses
    data is preferred. A formula found neither in gasses nor in materials project should be found in
    interpolation data, where the interpolations are done using materials project data.
    The use of interpolation is noted in each interpolation data entry because the interpolations give
    less precise dGf values than direct DFT calculations.
    Uses: used to determine dG_rxn at a given temperature
    """
    if CompTools(formula).n_els == 1:
        # elemental species have zero formation energy
        return {temp: 0 for temp in materials_project_data}
    try:
        formula = CompTools(formula).clean
    except:
        # if not cleanable, CompTools cannot be applied and dGrxn cannot be
        # calculated ***CHECK THIS***
        return None
    if formula in gas_data["300"]:
        return {temp: gas_data[temp][formula]["Ef"] for temp in gas_data}
    # first check gasses-- this data is preferred to MP if a formula is in both
    elif formula in materials_project_data["300"]:
        try:
            n_atoms = CompTools(formula).n_atoms
        except:
            return None
        if "C1O3" in formula or "CO3" in formula:
    #this correction accounts for DFT inaccuracy in calculating carbonates
            return {
                temp: ((materials_project_data[temp][formula]["Ef"]) * n_atoms + 0.830)
                / n_atoms
                for temp in materials_project_data
            }
        else:
            return {
                temp: materials_project_data[temp][formula]["Ef"]
                for temp in materials_project_data
            }
    elif not interpolation_data:
        # specifically allows for interpolation data modification (of no interpolation_data input)
        # because we want the interpolated formula to have zero energy for the purpose of
        # calculating dG_rxn for the interpolation 'reaction'-- see
        # modify_one_interpolation below
        return {temp: 0 for temp in materials_project_data}
    elif formula in interpolation_data:
        # note that the structure of the interpolation data database is different
        # from the other databases hence the different indexing
        return interpolation_data[formula]["dGf"]
    else:
        # if the formula is not found in any database, return None
        return None
    
def get_dGf_at_temp(
    formula, temp, materials_project_data, gas_data, original_interpolation_data
):
    """
    Args: formula string, temperature (string or float), databases
    Returns: dGf at the given temperature (interpolated as needed) or None if formula not found
    Function: This function is used to get the dGf of a formula at a specific temperature.
    For a temperature not in the database, the dGf is interpolated between the nearest temperatures.
    NOTE THE PROCEDURE FOR temperatures <300K and >2000K: for temperature below 300, 300 and 0 can be used
    if zero is a key. If not, 300 and 400 are used. For temperatures above 2000, 1900 and 2000 are used.
    Note that these interpolations are likely less accurate since the temperature may not lie between
    the two interpolation points, or may not be as close to either.
    Uses: used to determine dG_rxn at a given temperature from dGfs of all
    participating species
    """
    if type(temp) == str:
        temp = float(temp)
    dGf_dict = get_dGf_from_source(
        formula, materials_project_data, gas_data, original_interpolation_data
    )
    if dGf_dict == None:
        return None
    temps = [t for t in dGf_dict]
    if str(temp) in temps:
        return dGf_dict[str(temp)]
    # if temperature is in the database, return the dGf at that temperature
    # in that case, no temperature interpolation is required
    upper_T = math.ceil(temp / 100) * 100
    lower_T = upper_T - 100
    # finding the closes two temperatures in the database
    if temp >= 2000:
        lower_T = 1900
        upper_T = 2000
    # the maximum temperature is 2000, so for a temperature at or above 2000,
    # the interpolation is done between 1900 and 2000-- this makes the determination
    # slightly less accurate (especially for temperatures well above 2000)
    if temp <= 300:
        if "0" in temps:
            lower_T = 0
        else:
            lower_T = 300
        upper_T = 400
    # the minimum temperature is sometimes '0' and sometimes '300', but there are
    # no entries for temperatures in between, so for a temperature between 0 and 300,
    # if both keys are available, those two temperatures are used. If not, 300 and 400
    # are used. This makes the determination slightly less accurate (especially for very
    # low temperatures)
    temperatures_to_grab = [str(lower_T), str(upper_T)]
    dGfs = [dGf_dict[t] for t in temperatures_to_grab]
    if not dGfs:
        return None
    if dGfs[0] == None or dGfs[1] == None:
        return None
    slope = (dGfs[1] - dGfs[0]) / (
        float(temperatures_to_grab[1]) - float(temperatures_to_grab[0])
    )
    # interpolation-- difference in dGf divided by difference in temperature
    return dGfs[0] + slope * (temp - float(temperatures_to_grab[0]))
    # using point-slope form of a line to determine dGf at the given temperature

def get_chempot_correction(element, temp, env = 'air'):
    """
    Get the normalized correction term Δμ for chemical potential of a gas
    phase consisting of element at given temperature and pressure,
    referenced to that in the standard state (T_std = 298.15 K,
    T_std = 1 bar). The gas phase is limited to be one of O2, N2, Cl2,
    F2, H2. Calculation formula can be found in the documentation of
    Materials Project website.

    Args:
        element (string): The string representing the element.
        temp (float): The temperature of the gas phase.
        pres (float): The pressure of the gas phase.

    Returns:
        The correction of chemical potential in eV/atom of the gas
        phase at given temperature and pressure.
    """

    EV_TO_KJ_PER_MOL = 96.4853
    element = CompTools(element).clean
    # if not 'env':
    #     env = 'air'
    # if env == 'alumina':
    #     env = 'air'
    if element not in ['O', 'H', 'C1O2', 'H2O1','O2']:
        return 0
    if element == 'O2':
        element = 'O'
    std_temp = 298.15
    std_pres = 1E5
    ideal_gas_const = 8.3144598
    # Cp and S at standard state in J/(K.mol). Data from
    # https://janaf.nist.gov/tables/O-029.html
    # https://janaf.nist.gov/tables/N-023.html
    # https://janaf.nist.gov/tables/Cl-073.html
    # https://janaf.nist.gov/tables/F-054.html
    # https://janaf.nist.gov/tables/H-050.html
    Cp_dict = {'O': 29.376,
               'C1O2': 37.129,
               'H2O1': 33.22}
    S_dict = {'O': 205.147,
              'C1O2': 213.79,
              'H2O1': 194.10}
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
    
    Cp_std = Cp_dict[element]
    S_std = S_dict[element]
    pres = pp_dict[element]
    PV_correction = ideal_gas_const * temp * np.log(pres / std_pres)
    TS_correction = - Cp_std * (temp * np.log(temp) - std_temp * np.log(std_temp)) \
        + Cp_std * (temp - std_temp) * (1 + np.log(std_temp)) \
        - S_std * (temp - std_temp)

    dG = PV_correction + TS_correction

    # Convert to eV/molecule unit.
    dG /= 1000 * EV_TO_KJ_PER_MOL

    # Normalize by number of atoms in the gas molecule
    if element == 'H2O':
        dG /= 3
    if element == 'C1O2':
        dG /= 3
    # if element == 'NH3':
    #     dG /= 4
    if element == 'O':
        dG /= 2

    return dG

def get_amts_vars_target(reaction):
    """Args: reaction_data[i] for integer index i
    Returns: the reaction string associated with reactions of variable composition,
    which should include variable coefficients (or element amounts), or None if the
    reaction does not have variable composition
    This reaction string may be employed in the determination of entropy of disorder
    """
    target = reaction['target']['material_formula']
    vars = reaction['reaction']['thermo'][0]['amts_vars']
    if vars:
        return target
    else:
        return None

def parse_formula_for_amts_vars(amts_vars_target):
    # Pattern to find all element symbols possibly followed by numbers, variables, or algebraic expressions
    elements = [
    'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr',
    'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd',
    'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
    'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm',
    'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc', 'Lv', 'Ts', 'Og']

    pattern = '(' + '|'.join(sorted(elements, key=len, reverse=True)) + ')([0-9]*\.?[0-9]*[+\-]*[a-z0-9]*)'
    matches = re.findall(pattern, amts_vars_target)
    composition = {}

    for element, amount in matches:
        if amount == '':
            amount = 1  # Default coefficient is 1 if not specified
        elif re.match(r'^[\d.]+$', amount):  # Matches purely numeric coefficients
            amount = float(amount) if '.' in amount else int(amount)
            # If the amount is non-numeric, it remains as a string
            if type(amount) == str:
                amount = re.sub(r'(\d)([a-zA-Z])', r'\1*\2', amount)
            # Insert multiplication sign between a variable and a number
                amount = re.sub(r'([a-zA-Z])(\d)', r'\1*\2', amount)
        composition[element] = amount
    
    return composition

def get_target_sum(amts_vars,target_comp, unsubbed = False):
    for key, value in amts_vars.items():
        if unsubbed:
            exec(f"{key} = 0")
        else:
            exec(f"{key} = {value}")
            # print('key=', key)
    var_els = {}
    # print('target_comp', target_comp)
    for k,v in target_comp.items():
        if type(v) == str:
            try:
                v_new = float(eval(v))
                # print('v_new=', v_new)
                var_els[k] = v_new
            except:
                return None, None
    return var_els, sum(list(var_els.values()))

def get_var_target_comp(target_comp, amts_vars):
    var_target_comp = {var:[] for var in amts_vars}
    for var in amts_vars:
        for el in target_comp:
            if type(target_comp[el]) == str:
                if var in target_comp[el]:
                    var_target_comp[var].append(el)
    return var_target_comp


def get_amts_vars_correction(amts_vars,amts_vars_target,target):
    #Assume 1 site, skip cases that have no entropy 
    R = 8.617333262145E-5
    if not amts_vars:
        return 0
    # print(target)
    try:
        n_atoms = CompTools(target).n_atoms
    except:
        return None
    target_comp = parse_formula_for_amts_vars(amts_vars_target)
    # print('target_comp',target, target_comp,amts_vars_target)
    #dictionary of the form {element: coefficient} for the target
    sub_var, unsubbed_sum = get_target_sum(amts_vars,target_comp, unsubbed = True)
    # dictionary of the form {element: coefficient} evaluated with variables set to 0
    # and sum of coefficients (associated with variables) evaluated with variables set to 0
    var_els, var_sum = get_target_sum(amts_vars,target_comp)
    # same as above but with variables evaluated
    if var_sum != unsubbed_sum:
        return None
    if var_els:
        if all(amts_vars[var] == 0 for var in amts_vars):
            return 0
    var_target_comp = get_var_target_comp(target_comp, amts_vars)
    # print('var_target_comp',var_target_comp)
    if not var_els or not sub_var:
        # print('var els or sub var not found',amts_vars_target)
        return None
    # dictionary of the form {variable: [elements]} for the target-- all elements
    # that have coefficients associated with a given variable
    sub_sites = [el for el in var_els if var_els[el]<sub_var[el]]
    # print(amts_vars_target, sub_sites)
    if len(sub_sites) == 0:
        # print('no sub sites',amts_vars_target, amts_vars)
        return None
    # identify each element that has a coefficient that decreases with increase in
    # variable values
    site_dict = {el: 0 for el in var_els}
    # dictionary of the form {variable: site_count} for the target
    # here, I will consider site count to be the original number of atoms at the subbed site
    for el in var_els:
        if el in sub_sites:
            site_dict[el] = sub_var[el]
        else:
            sites = []
            vars = [var for var in var_target_comp if el in var_target_comp[var]]
            # print(amts_vars_target,'vars',vars)
            for var in vars:
                for el in var_target_comp[var]:
                    if el in sub_sites:
                        sites.append(sub_var[el]) 
            # print(amts_vars_target,'sites',sites)
            if len(sites) == 1:
                site_dict[el] = sites[0]
            else:
                # print('multiple sites found for substitution',amts_vars_target)
                return None
                
    #check which things sum to an integer 
    entropy = 0
    for v in var_els:
        if site_dict[v] == 0:
            return None
        entropy += -R*var_els[v]*np.log(var_els[v]/site_dict[v])/n_atoms
    return entropy


def get_ReactionEnergy(
    reaction_dict, temp,  materials_project_data, 
    gas_data, interpolation_data, env=None, 
    amts_vars = None, amts_vars_target = None, target = None,
    corrected = False
):
    """
    Args: reaction dictionary, temperature (string or float), databases
    Returns: ReactionEnergy object
    Uses: ReactionEnergy object can be used for various applications,
    like calculation dG_rxn and balancing a reaction
    """
    if not reaction_dict:
        return None
    if type(reaction_dict) != dict:
        return None
    if not temp:
        temp = 1073.15
    reactants = reaction_dict["reactants"]
    products = reaction_dict["products"]
    formulas = reactants + products
    dGf_dict = {str(formula): {"E": 0} for formula in formulas}
    # initialize dictionary of energies
    for formula in formulas:
        dGf = get_dGf_at_temp(
            formula, temp, materials_project_data, gas_data, interpolation_data
        )
        if dGf == None:
            return None
        if corrected:
            if formula in ['O2', 'H2O1', 'C1O2']:
                dGf_dict[formula]["E"] = dGf + get_chempot_correction(formula, temp,env)
            elif amts_vars:
                if formula in products:
                    amts_vars_correction = get_amts_vars_correction(amts_vars,amts_vars_target,target)
                    if amts_vars_correction == None:
                    # 'no correction found for amts_vars'
                        # print('no correction found for amts_vars')
                        return None
                    dGf_dict[formula]["E"] = dGf - temp*amts_vars_correction
                else:
                    dGf_dict[formula]["E"] = dGf
            else:
                dGf_dict[formula]["E"] = dGf
        else:
            dGf_dict[formula]["E"] = dGf
    re = ReactionEnergy(
        input_energies=dGf_dict,
        reactants=reactants,
        products=products,
        norm="atom",
        # dG_rxn is returned in ev/atom
        energy_key="E",
    )
    try:
        dG_rxn = re.dE_rxn
    except:
        return None
    return re

def get_balanced_reaction_coefficients(precursors, target):
    formulas = precursors + [target]
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



def get_dGrxn_at_T(
    reaction_dict, temp, materials_project_data, 
    gas_data, interpolation_data, corrected = False, 
    environment = None, amts_vars = None, 
    amts_vars_target = None, target = None
):
    """
    Args: reaction dictionary, temperature (string or float), databases
    Returns: dG_rxn at the given temperature (interpolated as needed)
    or None if formula not found
    Function: ReactionEnergy class is applied to reaction to
    extract class attribute dE_rxn
    """
    re = get_ReactionEnergy(
        reaction_dict, temp, materials_project_data, 
        gas_data, interpolation_data, corrected, environment, 
        amts_vars, amts_vars_target, target
)
    try:
        re.rxn_string
    except:
        print('no balanced reaction')
    try:
        dG_rxn = re.dE_rxn
        rxn = re.rxn_string
    except:
        print("cannot calculate dG_rxn from ReactionEnergy for", reaction_dict)
        return None, None
    return dG_rxn, rxn

# def get_energy_from_mp(elems,temperature):
#     mpr = MPRester(api_key)
#     # mp_entries = mpr.get_entries_in_chemsys(elems,additional_criteria={'energy_above_hull': (0,0.05)})
#     mp_entries = mpr.get_entries_in_chemsys(elems)
#     gcse =  GibbsComputedStructureEntry.from_entries(mp_entries, temperature)
#     # print([entry.as_dict()['formation_enthalpy_per_atom'] for entry in gcse])
#     return gcse

def get_gibbs_entry_set(elems,temperature):
    mpr = MPRester(api_key)
    mp_entries = mpr.get_entries_in_chemsys(elems)
    gibbs_entry_set = GibbsEntrySet.from_computed_entries(mp_entries, temperature,include_nist_data=True)
    gibbs_entry_set = gibbs_entry_set.filter_by_stability(e_above_hull=0.05)
    return gibbs_entry_set


def get_ComputedReaction(precursors, target, temperature, mp_data, gas_data, interpolation_data,query = False, gibbs_set = None):
    # print('cr target',target)
    coef_dict = get_balanced_reaction_coefficients(precursors, target)
    if not coef_dict:
        return None
    coef_dict = {key:coef_dict[key] for key in coef_dict if coef_dict[key]!=0}
    species = list(coef_dict.keys())
    # print('species',species)
    coefficients = list(coef_dict.values())
    comp_entries = []
    if not temperature:
        temperature = 1073.15
    for spec in species:
        # print(spec)
        n_atoms = CompTools(spec).n_atoms
        if not query:
            spec_energy = n_atoms*get_dGf_at_temp(spec,temperature, mp_data, gas_data, interpolation_data)
        else:
            elems = CompTools(spec).els
            if not gibbs_set:
                gibbs_set = get_gibbs_entry_set(elems,temperature)
            spec_entry = [entry for entry in gibbs_set if CompTools(entry.composition.reduced_formula).clean == CompTools(spec).clean]
            if not spec_entry:
                # print('no entry found for',spec)
                continue
            # spec_energy = n_atoms*spec_entry[0].formation_enthalpy_per_atom
        # print(spec,spec_energy)
        # comp = Composition(spec)
        # comp_entry = ComputedEntry(comp, spec_energy)
        comp_entries.append(spec_entry[0])
    cr = CR(comp_entries, coefficients,balanced = True)
    # print('computed reaction',cr)
    return cr

def get_chemsys_mp_targets_on_hull(mp_data, target, query = False):
    chemsys = CompTools(target).chemsys
    formulas = sorted(list(mp_data["300"].keys()))
    chemsys_formulas = []
    for formula in formulas:
        if CompTools(formula).chemsys == chemsys:
            if mp_data['300'][formula]['Ed'] < 0.05:
                chemsys_formulas.append(CompTools(formula).clean)
    if CompTools(target).clean not in chemsys_formulas:
        chemsys_formulas.append(target)
    return chemsys_formulas

def main():
    return

if __name__ == "__main__":
    main()