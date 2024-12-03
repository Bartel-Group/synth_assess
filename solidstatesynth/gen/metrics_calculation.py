import os
from pydmclab.utils.handy import read_json
from pydmclab.core.comp import CompTools

from rxn_network.core import Composition
from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.enumerators.basic import BasicEnumerator
from rxn_network.enumerators.basic import BasicOpenEnumerator
from rxn_network.reactions.hull import InterfaceReactionHull
from collections.abc import Iterable

from pydmclab.core.comp import CompTools

from pymatgen.entries.computed_entries import CompositionEnergyAdjustment
import math
from solidstatesynth.gen.entries import GibbsSet, FormulaChecker
import pickle

#Gas partial pressures in atm for different environments
PA_TO_ATM_CONV = 101300

GAS_PARTIAL_PRESSURES = {
    "air": {
        "O2": 21200 / PA_TO_ATM_CONV,
        "N2": 79033 / PA_TO_ATM_CONV,
        "CO2": 40.50 / PA_TO_ATM_CONV,#FIXED
        "H2O": 2300 / PA_TO_ATM_CONV,
        "H3N": 16 / PA_TO_ATM_CONV,
        "NO": 0.1 / PA_TO_ATM_CONV
    },
    "inert": {
        "O2": 0.1 / PA_TO_ATM_CONV,
        "N2": 0.1 / PA_TO_ATM_CONV,
        "CO2": 0.1 / PA_TO_ATM_CONV,
        "H2O": 0.1 / PA_TO_ATM_CONV,
        "H3N": 0.1 / PA_TO_ATM_CONV,
        "NO": 0.1 / PA_TO_ATM_CONV
    }
}


#boltzmann constant
KB = 8.6173303e-5


class EnumerateRxns():
    """
    Class for calculating competition metrics for single-step reactions in a single chemical system.
    Uses the reaction-networks package for calculating reaction energy, primary and secondary competition, and the gamma cost function.
    Can be utilized by imnputting select precursors, targets, and temperature, and will output the competition metrics for each reaction.
    INSERT CREDIT TO MATT MCDERMOTT + maybe others???
    
    """

    def __init__(
        self,
        els: Iterable[str],
        temperature: float = 300,
        with_theoretical: bool = True,
        stability_filter: float = 0.05,
        open = True,
        all_possible_precursors = None,
        remake: bool = True
    ):
        """
        Initialize the calculator with the specified elements and temperature
        Args:
            els (Iterable[str]): Formula Strings of all the elements in the chemical system.
            precursors (Iterable[str]): Optional way to specify the precursors to use in the reaction enumeration. Otherwise,
            the possible precursors will be identified from the elements and known text-mined precursors. Specified precursors
            will still be filtered by the elements in the chemical system.
        """

        self.els = list(set(els))
        if open:
            file_name = ''.join(els)+str(temperature)+ '_open_rxns.pkl'
        else:
            file_name = ''.join(els)+str(temperature) + '_rxns.pkl'
        if os.path.exists(os.path.join('../data/', file_name)) and remake == False:
            file = open('../data/' + file_name, 'rb')
            self.rxns = pickle.load(file)
        else:
            self.open = open
            self._temperature = temperature
            self.with_theoretical = with_theoretical
            self.stability_filter = stability_filter
            if not all_possible_precursors:
                all_possible_precursors = read_json(os.path.join('../data/', 'tm_precursors.json'))['data']
            self.all_possible_precursors = all_possible_precursors
            #Initialize the entries and reactions to None
            self.entries = None
            self.rxns = None
            self._build_calculator()
            with open('../data/' + file_name, 'wb') as f:
                pickle.dump(self.rxns, f)

        # numpy array should go deeper 


    def _get_entries(self):
        els = self.els
        temperature = self._temperature
        with_theoretical = self.with_theoretical
        stability_filter = self.stability_filter
        entry_set = GibbsSet(chemsys_els = els,temperature = temperature, with_theoretical=with_theoretical,stability_threshold=stability_filter).entry_set
        print('entry set',entry_set.entries)
        return entry_set

    def possible_precursors(self):
        els = self.els
        all_possible_precursors = self.all_possible_precursors
        possible_precursors = []
        for precursor in all_possible_precursors:
            if FormulaChecker(precursor, els).is_relevant:
                possible_precursors.append(precursor)
        return possible_precursors

    def _build_calculator(self):
        """
        Build the calculator by getting the entries and enumerating the reactions
        """
        self.entries = self._get_entries()
        entries = self.entries
        precursors = self.possible_precursors()

        print('entries obtained')
        #Use the BasicEnumerator and BasicOpenEnumerator to enumerate all the reactions from the precursors
        self.rxns = BasicEnumerator(precursors = precursors, exclusive_precursors=True).enumerate(entries)

        print('rxns enumerated')
        if self.open:
            self.rxns = self.rxns.add_rxn_set(BasicOpenEnumerator(open_phases=["O2"], precursors = precursors, exclusive_precursors=True).enumerate(entries))
        #Filter out duplicate reactions
        self.rxns = self.rxns.filter_duplicates()
        print('rxn duplicates filtered')
        return

    


class TargetRxns():
    """
    class for calculating competition
    """
    def __init__(self, target, reactions, temperature = 300, environment = "air", open = True):
        self.target = target
        self.reactions = reactions
        self.temperature = temperature
        self.environment = environment
        self.open = open


    def target_rxns(
            self,
            reactions = None   
        ) -> list[ComputedReaction]:
        """
        Helper method to find all the reactions that make the specified target without any byproducts (besides balancing gases)
        Args:
            target (str): String of the target to find the reactions for
            reactions (ReactionSet): ReactionSet of all the reactions to search through at the reaction temperature
        Returns:
            list[ComputedReaction]: List of all the reactions that make the target without any byproducts
        """
        if not reactions:
            reactions = self.reactions
        else:
            reactions = reactions
        target = self.target
        
        target_rxns = []

        #The products that are allowed in the target reactions - the target and any balancing gases
        allowed_products = {Composition(target)} | set([Composition(i) for i in GAS_PARTIAL_PRESSURES["air"].keys()])

        #Check each reaction to see if its products are within the allowed products
        #The precursors should already be the only reactants in the reaction - besides O2 in some cases
        #These are already filtered by the correct atmosphere - if it is open to air, O2 can be a reactant

        for rxn in reactions:
            if CompTools(target).clean in [CompTools(i).clean for i in rxn.products]:
                if set([i.reduced_composition for i in rxn.products]) <= allowed_products and (len(rxn.reactants) == 2 or len(set(rxn.reactants) - {Composition("O2")}) == 2):
                    target_rxns.append(rxn)
        print('target rxns found')
        print(target_rxns)

        return ReactionSet.from_rxns(target_rxns)
    

    def environment_correction(
            self,
            formula: str,
        ) -> CompositionEnergyAdjustment:
        """
        Get the environment correction for a gas at the specified temperature
        Different environments will have different partial pressures for gasses, which will affect the correction.
        Args:
            formula (str): Formula of the gas to get the environment correction for
            temperature (float): Temperature in Kelvin to get the environment correction for
            environment (str): Environment to correct the entry for, either "air" for air environment or "inert" for excess inert environment
        Returns:
            CompositionEnergyAdjustment: Environment correction
        """
        temperature = self.temperature
        environment = self.environment
        open = self.open
        if environment not in GAS_PARTIAL_PRESSURES:
            return None

        n_atoms = CompTools(formula).n_atoms
        #Calculate the environment correction for the gas on a per atom basis
        if not open:
            adjustment_per_atom = 0
        else:
            adjustment_per_atom = KB * math.log(GAS_PARTIAL_PRESSURES[environment][formula])*float(temperature)/n_atoms
        #Set a useful name for the correction - This is important because the correction does not automatically change with temperature
        name = f"{formula} {environment} correction @ {temperature}K"

        return CompositionEnergyAdjustment(adj_per_atom=adjustment_per_atom, n_atoms=n_atoms, name=name)


    def rxns_at_temp_env(
        self,
    ) -> ReactionSet:
        """
        Get the reactions at the specified temperature and environment. This requires adding a correction for the partial pressure of gases in the environment.
        For an air environment, O2 is a viable reactant, but not for an excess inert atmosphere.
        Args:
            temp (float): Temperature in Kelvin to get the reactions at
            environment (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            ReactionSet: Reactions at the specified temperature
        """
        temperature = self.temperature
        environment = self.environment
        open = self.open
        rxns = self.reactions
        #Chemical potential of O at specified temperature and evironment partial pressure
        mu = KB * (temperature) * math.log(GAS_PARTIAL_PRESSURES[environment]["O2"])

        #For inert environment, remove all reactions that contain O2 as a reactant
        #Determine the allowed gas products to add the correction to
        gas_comps = [i.composition for i in rxns.entries if i.reduced_formula in set(GAS_PARTIAL_PRESSURES["air"].keys())-{"O2"}]
        #Get the environment correction for each gas
        gas_dict = {i: self.environment_correction(i.reduced_formula) for i in gas_comps}
        
        # #Add the environment correction to the gas products
        for i in rxns.entries:
            if Composition(i.composition) in gas_comps:
                i.energy_adjustments.append(gas_dict[Composition(i.composition)])
                print('correction',i.composition, i.energy_adjustments)

        #Finally set O as an open element with the calculated chemical potential from partial pressure
        if not open:
            return rxns
        else:
            return rxns.set_chempot(open_el="O", chempot=mu)

    def metrics_at_temp_env(
        self,
        rxns = None
    ) -> list[dict[str, ComputedReaction | float]]:
        """
        Calculate the competition metrics for all the reactions at the specified temperature and environment
        Args:
            rxns (ReactionSet) : optional argument
        Returns:
            metrics (list[dict[str, ComputedReaction | float]]): List of dictionaries with the reaction, energy, primary competition, secondary competition, and gamma cost function for each reaction
        """

        environment = self.environment

        #Get the reactions at the correct temperature
        if not rxns:
            print('generating')
            rxns = self.rxns_at_temp_env()
        target_rxns_at_temp = self.target_rxns(rxns)
        

        #Initialize the metrics data
        metrics = []

        #Loop through all the targets to calculate metrics for all the precursors combinations to make each target.
        #Find the reactions that make the target without any byproducts - will calculate metrics for each reaction
        for rxn in target_rxns_at_temp:
            reactants = rxn.reactants
            if len(reactants) == 1:
                continue
            #Get the precursors for the reaction - besides O2 if there are 2 solid precursors
            rxn_prec = set([i.reduced_composition for i in rxn.reactants])
            if len(rxn_prec) == 3 and Composition("O2") in rxn_prec:
                rxn_prec.remove(Composition("O2"))
            print(rxn_prec)
            #Get the filtered reactions for the InterfaceReactionHull that only contain the precursors in this reaction
            #In air, O2 is allowed as a reactant
            #In excess inert, O2 is not allowed as a reactant
            if environment == "air":
                filtered_rxns = list(rxns.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]+["O2"]))
            else:
                filtered_rxns = list(rxns.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]))

            filtered_rxns = [i for i in filtered_rxns if len(i.reactants) > 1]
            #Build the InterfaceReactionHull from the precursors and filtered reactions

            irh = InterfaceReactionHull(*rxn_prec, filtered_rxns)

            # irh.plot().show()
            rxn_data = {}
            rxn_data["rxn"] = str(rxn)
            rxn_data["energy"] = rxn.energy_per_atom
            rxn_data["c1"] = irh.get_primary_competition(rxn)
            rxn_data["c2"] = irh.get_secondary_competition(rxn)
            rxn_data["gamma"] = 0.1 * rxn_data["energy"] + 0.45 * rxn_data["c1"] + 0.45 * rxn_data["c2"]
            
            # #Useful for debugging:
            # rxn_data["competing_rxns"] = [str(i)+" dG="+str(i.energy_per_atom) for i in filtered_rxns]

            metrics.append(rxn_data)

        return metrics

def get_metrics(target, temperature, open = True):
    r = EnumerateRxns(els = CompTools(target).els, temperature = temperature).rxns
    t = TargetRxns(target = target, reactions = r, temperature = temperature, open = open)
    return t.metrics_at_temp_env()   

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


def get_competitions(nist, with_nist = False):
    DATA_DIR = "/Users/schle759/Documents"
    fjson = os.path.join(DATA_DIR, "mcdermott.json")
    competition = read_json(fjson)
    # competitions = [entry for entry in competition if 'open' in entry]
    new_entries = []
    for entry in competition:
        rxn_dict = get_reaction_dict_from_string(entry['reaction'])
        temp = float(entry['temp_degC']) + 273.15
        if temp > 2000:
            temp = 2000
        if 'open' not in entry:
            open = False
        else:
            open = True
        entry_new = {'precursors': rxn_dict['reactants'],  
        'temp': temp, 'energy': entry['energy'],
        'c1': entry['c1'], 'c2': entry['c2'], 'open': open}
        target = [formula for formula in rxn_dict['products'] if CompTools(formula).clean not in ['O1', 'H2O1','C1O2']][0]
        entry_new['target'] = target
        if not with_nist:
            if not_nist_rxn(entry_new, nist):
                new_entries.append(entry_new)
        else:
            new_entries.append(entry_new)
    return new_entries

def nist_cmpds():
    DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"
    f = os.path.join(DATADIR, "nist_compounds.json")
    cmpds = read_json(f)
    nist_cmpds = [CompTools(entry).clean for entry in cmpds]
    return nist_cmpds

def not_nist_rxn(competition_entry, nist):
    # rxns = []
    species = competition_entry['precursors'] + [competition_entry['target']]
    if all([CompTools(i).clean not in nist for i in species]):
        return True
    return False
