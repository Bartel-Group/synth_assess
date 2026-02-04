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
from synth_assess.selectivity.entries import GibbsSet, FormulaChecker
from synth_assess.data.load import mp_data



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
"""
This code draws on the work of McDermott, M. J., Dwaraknath, S. S., and Persson, K. A.  https://doi.org/10.1038/s41467-021-23339-x
"""


class PrecursorSet():
    def __init__(
                 self, 
                 els: Iterable[str], 
                 solids_data: dict = None, 
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
            If not user-specified, mp data is used for solids_data (from which to generate possible reactions)
            restrict_to_tm (bool): Whether to restrict the precursors to only those listed in the text-mined dataset
            with_theoretical (bool): Whether to include theoretical MP entries in the precursor set
            filter_data (bool): Whether to filter the data for the specified elements and/or based on tm/theoretical args
            Note that, implicitly, if restrict_to_tm and with_theoretical are both False and filter_data = True,  the precursor set will include only 
            ICSD materials from MP. If filter_data = False, the list of precursors will include all materials in the solids data. 
        """
        self.els = els
        if solids_data:
            self.solids_data = solids_data
        else:
            self.solids_data = mp_data()
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
    Can be utilized by inputting select precursors, targets, and temperature, and will output the competition metrics for each reaction.
    Warning: the code to generate reaction networks will be parallelized over all available processors.
    
    """

    def __init__(
        self,
        els: Iterable[str],
        solids_data: dict = None,
        temperature: float = 300,
        gibbs_kwargs: dict = {},
        prec_kwargs: dict = {},
        gen_data = None,
        gen_formula = None,
    ):
        """
        Initialize the calculator with the specified elements and temperature
        Args:
            els (Iterable[str]): Formula Strings of all the elements in the chemical system.
            if not user-specified, MP data is used for solids_data.
            gibbs kwargs and prec kwargs are a dictionary of relevant arguments from GibbsSet and PrecursorSet respectively
        """

        self.els = list(set(els))
        self.temperature = temperature
        if solids_data:
            self.solids_data = solids_data
        else:
            self.solids_data = mp_data()
        self.temperature = temperature
        self.prec_kwargs = prec_kwargs
        self.gibbs_kwargs = gibbs_kwargs
        self.gen_data = gen_data
        self.gen_formula = gen_formula
        #Initialize the entries and reactions to None
        self.entries = None
        self.rxns = None
        self._build_calculator()



    def _get_entries(self):
        """refer to solidstatesynth.gen.entries GibbsSet class for kwargs (additional control to set)"""
        els = self.els
        kwargs = self.gibbs_kwargs
        solids_data = self.solids_data
        temperature = self.temperature
        entry_set = GibbsSet(chemsys_els = els,temperature = temperature,  
                             solids_data = solids_data, gen_data = self.gen_data,
                             gen_formula = self.gen_formula, **kwargs).entry_set
        return entry_set


    def _build_calculator(self):
        """
        Build the calculator by getting the entries and enumerating the reactions.
        Executing the calculator generates the set of rxns (callable by EnumerateRxns(...).rxns)
        """
        self.entries = self._get_entries()
        entries = self.entries
        kwargs = self.prec_kwargs
        precursors = PrecursorSet(els = self.els, solids_data=self.solids_data, **kwargs).precursors
        print('entries obtained')
        self.rxns = BasicEnumerator(precursors = precursors, exclusive_precursors=True).enumerate(entries)
        print('rxns enumerated')
        self.rxns = self.rxns.filter_duplicates()
        print('rxn duplicates filtered')
        return



class TempEnvCorrections():
    """
    Applying environment/temperature corrections based on variable gas partial pressure and variable chemical potential
    with respect to O. Note that these corrections are distinct from the SISSO temperature correction which is applied
    using the Gibbs entry object.
    """
    def __init__(self,  
                 temperature: float = 300, 
                 environment: str = "air",
                 ):
        self.temperature = temperature
        self.environment = environment

        
    def environment_correction(
            self,
            formula: str,
        ) -> CompositionEnergyAdjustment:
        """
        Get the environment correction for a gas at the specified temperature
        This will be zero for all cases of a "closed" reaction. This function is included only for the utility
        of the user if an open reaction computation is desired.

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
        if environment not in GAS_PARTIAL_PRESSURES:
            return None

        n_atoms = CompTools(formula).n_atoms
        adjustment_per_atom = 0
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
        gasses = set(GAS_PARTIAL_PRESSURES[environment].keys())-{"O2"}
        # #Add the environment correction to the gas products
        for i in rxns.entries:
            formula = i.reduced_formula
            if formula in gasses:
                i.energy_adjustments.append(self.environment_correction(formula))

        return rxns

class CustomClass(ReactionSet):
    """
    A helper subclass to modify McDermott's ReactionSet for our purposes:
    in this class, we redefine the "from_rxns" function in order to constrain entry
    IDs to facilitate matching entries at new temperatures to pre-generated reactions
    """
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

        # Always use entry_id, ignoring unique_id (the main distinction between this code
        # and the original code)
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

        # this reaction set allows for hard-coded "entry_ids"
        return rxn_set

class RxnsAtNewTempEnv():
    def __init__(self,
                 reaction_set: ReactionSet,
                 els: Iterable[str],
                 solids_data: dict = None ,
                 new_temperature: float = 300,
                 environment: str = "air",
                 gen_data = None,
                 gen_formula = None,
                 original_temperature = 300
                ):
        
        if solids_data:
            self.solids_data = solids_data
        else:
            self.solids_data = mp_data()
        self.temperature = new_temperature
        self.environment = environment
        self.corr = TempEnvCorrections(temperature = new_temperature,
                                       environment=environment)
        self.entries = reaction_set.entries
        self.reactions = reaction_set
        self.gen_data = gen_data
        self.gen_formula = gen_formula
        self.els = els
        # The default original temperature is 300 K: generally entries will first
        # be generated at 300 and then modified. However, the modification is feasible
        # at other temperatures as well and the function will be modified accordingly
        # with the designated original temperature
        self.original_temperature = original_temperature

    def reaction_entry_ids(self):
        """
        Returns a dictionary of entry ids associated with each formula and its entry:
        this function is used to match previously generated reactions to Gibbs entries
        constructed at a new temperature. The new-temperature entries will have an enforced
        entry id designated by the entry ids associated with the reaction set entries
        {formula (cleaned): entry_id from generated reaction entry (str)
        for formula in reactions (each formula corresponding to a distinct entry)}

        """
        reaction_entries = self.reactions.entries 
        entry_id_dict = {CompTools(entry.composition.reduced_formula).clean: 
                         entry.entry_id for entry in reaction_entries}
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
        entries = GibbsSet(chemsys_els= self.els, solids_data=solids_data, 
                            temperature=self.temperature, entry_id_dict= entry_id_dict,
                            gen_data = self.gen_data, gen_formula = self.gen_formula).entries
        for entry in entries:
            entry_id = entry.entry_id
            if 'Experimental' in entry_id:
                # ExperimentalReferenceEntry type entries do not allow for a 
                # specified entry id and explicitly include temperature. Here,
                # we modify the temperature so that the entry id matches that of
                # the original entries, generated at 300 K. 
                entry.entry_id = entry_id.split('_')[0]+'_' + str(self.original_temperature)

        reactions = list(self.reactions)
        open_elem = None
        chempot = 0.0

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
                 environment: str = "air"):
        
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
            #Get the filtered reactions for the InterfaceReactionHull that only contain the precursors in this reaction
            #In air, O2 is allowed as a reactant
            #In excess inert, O2 is not allowed as a reactant
            if environment == "air":
                filtered_rxns = list(rxns.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]+["O2"]))
                # print('f',filtered_rxns)
            else:
                filtered_rxns = list(rxns.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]))

            filtered_rxns = [i for i in filtered_rxns if len(i.reactants) > 1]
            #Build the InterfaceReactionHull from the precursors and filtered reactions

            irh = InterfaceReactionHull(*rxn_prec, filtered_rxns)
            rxn_data = {}
            rxn_data["rxn"] = str(rxn)
            rxn_data["energy"] = rxn.energy_per_atom
            rxn_data["c1"] = irh.get_primary_competition(rxn)
            rxn_data["c2"] = irh.get_secondary_competition(rxn)
            # gamma here employs are new weights. weights can be user specified.
            rxn_data["gamma"] = 0.16 * rxn_data["energy"] + 0.56 * rxn_data["c1"] + 0.28 * rxn_data["c2"]
            
            metrics.append(rxn_data)

        return metrics
    

    
class GammaFromTarget():
    def __init__(self,
                 target: str = None,
                 temperature: float = 1073,
                 solids_data: dict = None):
        """
        This class implements the entire pipeline from input target/temperature to output metrics.
        Note that by default, the solids data used is MP data.
        Warning: the code to generate reaction networks will be parallelized over all available processors.
        """
        self.target = CompTools(target).clean
        self.temperature = temperature
        if solids_data:
            self.solids_data = solids_data
        else:
            self.solids_data = mp_data()
        
    
    def get_metrics(self, gen_data = None, is_gen = False, restrict_to_tm = False):
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
        target = self. target
        temperature = self.temperature
        solids_data = self.solids_data
        if is_gen:
            gen_formula = target
        else:
            gen_formula = None
        r = EnumerateRxns(els = CompTools(target).els,temperature = 300, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula, prec_kwargs={'restrict_to_tm_precursors':restrict_to_tm}).rxns
        r = RxnsAtNewTempEnv(reaction_set = r, els = CompTools(target).els, new_temperature = temperature, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula).corrected_reactions_at_temp()
        t = AnalyzeReactionSet(target = target, reactions = r, temperature = temperature)
        return t.metrics_at_temp_env()
    
    def opt_rxn(self, gen_data = None, is_gen = False,):
        """
        Identify and return the lowest-gamma (most selective) reaction entry. 
        """
        all_rxn_metrics = self.get_metrics(gen_data=gen_data, is_gen = is_gen)
        min_entry = min(all_rxn_metrics, key=lambda d: d['gamma'])
        return min_entry