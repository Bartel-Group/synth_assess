import os
import math
import numpy as np
from collections.abc import Iterable
from pymatgen.entries.computed_entries import CompositionEnergyAdjustment
from pymatgen.core.composition import Composition
from pydmclab.utils.handy import read_json
from pydmclab.core.comp import CompTools
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.enumerators.basic import BasicEnumerator, BasicOpenEnumerator
from rxn_network.reactions.hull import InterfaceReactionHull
from rxn_network.reactions.reaction_set import ReactionSet
from solidstatesynth.gen.entries import GibbsSet, FormulaChecker


#Gas partial pressures in atm for different environments
PA_TO_ATM_CONV = 101300

GAS_PARTIAL_PRESSURES = {
    "air": {
        "O2": 21200 / PA_TO_ATM_CONV,
        "N1": 79033 / PA_TO_ATM_CONV,
        "C1O2": 40.50 / PA_TO_ATM_CONV,#FIXED
        "H2O1": 2300 / PA_TO_ATM_CONV,
        "H3N1": 16 / PA_TO_ATM_CONV,
        "N1O1": 0.1 / PA_TO_ATM_CONV
    },
    "inert": {
        "O2": 0.1 / PA_TO_ATM_CONV,
        "N1": 0.1 / PA_TO_ATM_CONV,
        "C1O2": 0.1 / PA_TO_ATM_CONV,
        "H2O1": 0.1 / PA_TO_ATM_CONV,
        "H3N1": 0.1 / PA_TO_ATM_CONV,
        "N1O1": 0.1 / PA_TO_ATM_CONV
    }
}


#boltzmann constant
KB = 8.6173303e-5
DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"
DATADIR_enumerate = "/Volumes/cems_bartel/projects/negative-examples/data/rxn_networks"

class PrecursorSet():
    def __init__(
                 self, 
                 els: Iterable[str], 
                 solids_data: dict,
                 filter_precursors: bool = True,
                 precursor_stability_filter: float = 0.05, 
                 restrict_to_tm_precursors: bool = True, 
                 with_theoretical_precursors: bool = False,
                 ):
        """
        Initialize the precursor set with the specified elements
        Args:
            els (Iterable[str]): Formula Strings of all the elements in the chemical system.
            solids_data dict: Dictionary of the solids data to use for the precursors, of the form {formula:{data for formula}}
            -- for filtering, data must include keys "theoretical" (bool), "textmined_precursor" (bool) and "energy_above_hull" (float)
            restrict_to_tm (bool): Whether to restrict the precursors to only those listed in the text-mined dataset
            with_theoretical (bool): Whether to include theoretical MP entries in the precursor set
            filter_data (bool): Whether to filter the data for the specified elements and/or based on tm/theoretical args
            Note that, implicitly, if restrict_to_tm and with_theoretical are both False and filter_data = True,  the precursor set will include only 
            ICSD materials from MP. If filter_data = False, the list of precursors will include all materials in the solids data. 
        """
        self.els = els
        self.solids_data = solids_data
        self.filter_data = filter_precursors
        self.stability_filter = precursor_stability_filter
        self.restrict_to_tm = restrict_to_tm_precursors
        self.with_theoretical = with_theoretical_precursors
        
        
    def relevant_precursor(self, formula, entry):
        """
        each entry corresponds to an entry in solids data
        
        """
        if self.restrict_to_tm and not entry['tm_precursor']:
            return False
        if not self.with_theoretical and entry['theoretical']:
            return False
        if not FormulaChecker(formula, self.els).is_relevant:
            return False
        return True
    
    @property
    def precursors(self):
        if not self.filter_data:
            # to use the whole data set as possible precursors
            return list(self.solids_data.keys())
        precursors = []
        for formula, entry in self.solids_data.items():
            if self.relevant_precursor(formula, entry):
                precursors.append(formula)
        return precursors
    

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
        solids_data: dict,
        temperature: float = 300,
        open = True,
        gibbs_kwargs: dict = {},
        prec_kwargs: dict = {},
        remake: bool = True,
    ):
        """
        Initialize the calculator with the specified elements and temperature
        Args:
            els (Iterable[str]): Formula Strings of all the elements in the chemical system.
            gibbs kwargs and prec kwargs are a dictionary of relevant arguments from GibbsSet and PrecursorSet respectively
        """

        self.els = list(set(els))
        self.temperature = temperature
        # if open:
        #     file_name = ''.join(els)+str(temperature)+ '_open_rxns.pkl'
        # else:
        #     file_name = ''.join(els)+str(temperature) + '_rxns.pkl'
        # if os.path.exists(os.path.join(DATADIR_enumerate, file_name)) and remake == False:
        #     file = open(DATADIR_enumerate + file_name, 'rb')
        #     self.rxns = pickle.load(file)
        # else:
        self.solids_data = solids_data
        self.temperature = temperature
        self.open = open
        self.prec_kwargs = prec_kwargs
        self.gibbs_kwargs = gibbs_kwargs
        #Initialize the entries and reactions to None
        self.entries = None
        self.rxns = None
        self._build_calculator()
            # with open(DATADIR_enumerate + file_name, 'wb') as f:
            #     pickle.dump(self.rxns, f)
            
        # numpy array should go deeper 


    def _get_entries(self):
        els = self.els
        kwargs = self.gibbs_kwargs
        solids_data = self.solids_data
        temperature = self.temperature
        entry_set = GibbsSet(chemsys_els = els,temperature = temperature,  
                             solids_data = solids_data, **kwargs).entry_set
        return entry_set


    def _build_calculator(self):
        """
        Build the calculator by getting the entries and enumerating the reactions
        """
        self.entries = self._get_entries()
        entries = self.entries
        kwargs = self.prec_kwargs
        precursors = PrecursorSet(els = self.els, solids_data=self.solids_data , **kwargs).precursors

        print('entries obtained')
        #Use the BasicEnumerator and BasicOpenEnumerator to enumerate all the reactions from the precursors
        self.rxns = BasicEnumerator(precursors = precursors, exclusive_precursors=True).enumerate(entries)

        print('rxns enumerated')
        if self.open:
            self.rxns = self.rxns.add_rxn_set(BasicOpenEnumerator(open_phases=["O2"], 
                                                                  precursors = precursors, 
                                                                  exclusive_precursors=True).enumerate(entries))
        #Filter out duplicate reactions
        self.rxns = self.rxns.filter_duplicates()
        print('rxn duplicates filtered')
        return

#start with the same solids data 


class TempEnvCorrections():
    """
    Applying environment/temperature corrections based on variable gas partial pressure and variable chemical potential
    with respect to O. Note that these corrections are distinct from the SISSO temperature correction which is applied
    using the Gibbs entry object.
    """
    def __init__(self,  
                 temperature: float = 300, 
                 environment: str = "air",
                 open = True
                 ):
        self.temperature = temperature
        self.environment = environment
        self.open = open
        
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
        elif formula not in GAS_PARTIAL_PRESSURES[environment]:
            adjustment_per_atom = 0
        else:
            if CompTools(formula).clean == ['O1']:
                formula = 'O2'
            adjustment_per_atom = KB * math.log(GAS_PARTIAL_PRESSURES[environment][formula])*float(temperature)/n_atoms
        #Set a useful name for the correction - This is important because the correction does not automatically change with temperature
        name = f"{formula} {environment} correction @ {temperature}K"

        return CompositionEnergyAdjustment(adj_per_atom=adjustment_per_atom, n_atoms=n_atoms, name=name)

    @property
    def mu(self):
        """
        Get the chemical potential of O at the specified temperature and environment partial pressure
        Args:
            temp (float): Temperature in Kelvin to get the chemical potential at
            environment (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            float: Chemical potential of O at the specified temperature
        """
        temperature = self.temperature
        environment = self.environment
        if environment not in GAS_PARTIAL_PRESSURES:
            return None
        mu = KB * (temperature) * math.log(GAS_PARTIAL_PRESSURES[environment]["O2"])
        return mu
    
    def rxns_with_temp_env_correction(
        self,
        rxns: ReactionSet
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
        environment = self.environment
        open = self.open
        gasses = set(GAS_PARTIAL_PRESSURES[environment].keys())-{"O2"}
        #Chemical potential of O at specified temperature and evironment partial pressure

        #For inert environment, remove all reactions that contain O2 as a reactant
        #Determine the allowed gas products to add the correction to
        #Get the environment correction for each gas
        # gas_dict = {i: corrections.environment_correction(i.reduced_formula) for i in gas_comps}
        
        # #Add the environment correction to the gas products
        for i in rxns.entries:
            formula = i.reduced_formula
            if formula in gasses:
                i.energy_adjustments.append(self.environment_correction(formula))

        #Finally set O as an open element with the calculated chemical potential from partial pressure
        if not open:
            return rxns
        else:
            return rxns.set_chempot(open_el="O", chempot=self.mu)
        # what happens if i set mu to 'None' instead of zero? 1111\

class CustomClass(ReactionSet):
    @classmethod
    def from_rxns(
        cls,
        rxns: [ComputedReaction],
        entries: [GibbsSet],
        open_elem: str | None = None,
        chempot: float = 0.0,
        filter_duplicates: bool = False,
    ) -> ReactionSet:
        
        entries = sorted(set(entries), key=lambda r: r.composition)

        # Always use entry_id, ignoring unique_id
        all_entry_indices = {entry.entry_id: idx for idx, entry in enumerate(entries)}
        indices, coeffs, data = {}, {}, {}  # type: ignore

        for rxn in rxns:
            size = len(rxn.entries)

            rxn_indices = []
            for e in rxn.entries:
                rxn_indices.append(all_entry_indices[e.entry_id])
        
            if size not in indices:
                indices[size] = []
                coeffs[size] = []
                data[size] = []

            indices[size].append(rxn_indices)
            coeffs[size].append(rxn.coefficients)
            data[size].append(rxn.data)

        for size in indices:
            indices[size] = np.array(indices[size])
            coeffs[size] = np.array(coeffs[size])
            data[size] = np.array(data[size])

        all_open_elems: set[Element] = set()
        all_chempots: set[float] = set()

        if all(r.__class__.__name__ == "OpenComputedReaction" for r in rxns) and not open_elem:
            for r in rxns:
                all_open_elems.update(r.chempots.keys())
                all_chempots.update(r.chempots.values())

            if len(all_chempots) == 1 and len(all_open_elems) == 1:
                chempot = all_chempots.pop()
                open_elem = all_open_elems.pop()

        rxn_set = cls(
            entries=entries,
            indices=indices,
            coeffs=coeffs,
            open_elem=open_elem,
            chempot=chempot,
            all_data=data,
        )

        if filter_duplicates:
            rxn_set = rxn_set.filter_duplicates()

        return rxn_set
        # Optionally include logging or checks for debugging

class RxnsAtNewTempEnv():
    def __init__(self,
                 reaction_set: ReactionSet,
                 els: Iterable[str],
                 solids_data: dict ,
                 new_temperature: float = 300,
                 environment: str = "air",
                 open = True, 
                ):
        self.temperature = new_temperature
        self.environment = environment
        self.corr = TempEnvCorrections(temperature = new_temperature,
                                       environment=environment, 
                                       open = open)
        self.entries = reaction_set.entries
        self.reactions = reaction_set
        self.solids_data = solids_data
        self.els = els
        self.open = open

    def reaction_entry_ids(self):
        reaction_entries = self.reactions.entries 
        entry_id_dict = {CompTools(entry.composition.reduced_formula).clean: 
                         entry.entry_id for entry in reaction_entries}
        # for key in entry_id_dict:
        #     current_entry_id = entry_id_dict[key]
        #     if 'Experimental' in current_entry_id:
        #         id_no_temp = current_entry_id.split('_')[0]
        #         id_new = id_no_temp + f"_{self.temperature}"
        #         # id_new = id_no_temp + "_300"

        #         entry_id_dict[key] = id_new
        #         # print('id_new',id_new)
        # print(entry_id_dict)
        return entry_id_dict
    

    def reactions_at_temp(self):
        """
        Get the reactions at the specified temperature
        Args:
            temp (float): Temperature in Kelvin to get the reactions at
        Returns:
            ReactionSet: Reactions at the specified temperature
        """
        solids_data = self.solids_data
        entry_id_dict = self.reaction_entry_ids()
        # for key in entry_id_dict:
        #     if '_300' not in entry_id_dict[key]:
        #         entry_id_dict[key] = entry_id_dict[key] + '_300'
        print(entry_id_dict)
        # entries regenerated at new temperatures
        entries = GibbsSet(chemsys_els= self.els, solids_data=solids_data, temperature=self.temperature, entry_id_dict= entry_id_dict).entries
        for entry in entries:
            entry_id = entry.entry_id
            if 'Experimental' in entry_id:
                entry.entry_id = entry_id.split('_')[0]+'_300'
        # print(entries)
        print([entry.entry_id for entry in entries])
        print([entry.entry_id for entry in self.reactions.entries])
        # print(self.reactions.entries)
        reactions = list(self.reactions)
        if self.open:
            open_elem = 'O'
            chempot = self.corr.mu
        else:
            open_elem = None
            chempot = 0.0

        # print(list(reactions))
        rxns_at_temp = CustomClass.from_rxns(rxns = reactions, entries = entries, 
                                             open_elem= open_elem, chempot= chempot)
        return rxns_at_temp
    
    def corrected_reactions_at_temp(self):
        rxns = self.reactions_at_temp()
        corr = self.corr
        return corr.rxns_with_temp_env_correction(rxns)

    
class AnalyzeReactionSet():
    def __init__(self,
                 reactions: ReactionSet,
                 target: str = None,
                 temperature: float = 300,
                 environment: str = "air",):
        
        self.reactions = reactions
        self.temperature = temperature
        self.environment = environment
        self.target = target
    
    def target_rxns(
            self):
        """
        Get the reactions that make the target without any byproducts
        Args:
            target (str): Formula of the target compound
        Returns:
            list[ComputedReaction]: Reactions that make the target without any byproducts
        """ 
        target = self.target
        environment = self.environment
        rxns = self.reactions
        rxns_with_target = []
        allowed_products = {Composition(target)} | set([Composition(i) 
                                for i in GAS_PARTIAL_PRESSURES[environment].keys()])
        for rxn in rxns:
            if CompTools(target).clean in [CompTools(i).clean for i in rxn.products]:
                if set([i.reduced_composition for i in rxn.products]) <= allowed_products and (len(rxn.reactants) == 2 or len(set(rxn.reactants) - {Composition("O2")}) == 2):
                    rxns_with_target.append(rxn)
        print('target rxns found')

        return ReactionSet.from_rxns(rxns_with_target)


    def metrics_at_temp_env(
        self,
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

        rxns = self.reactions
        target_rxns_at_temp = self.target_rxns()
        

        #Initialize the metrics data
        metrics = []

        #Loop through all the targets to calculate metrics for all the precursors combinations to make each target.
        #Find the reactions that make the target without any byproducts - will calculate metrics for each reaction
        #correct entries based on atmosphere 
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
