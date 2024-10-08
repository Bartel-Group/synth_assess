import os
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools

from rxn_network.core import Composition
from pymatgen.core import Composition as PymatgenComposition
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.reactions.reaction_set import ReactionSet
from rxn_network.reactions.computed import ComputedReaction
from rxn_network.reactions.open import OpenComputedReaction

from pymatgen.entries.computed_entries import ComputedStructureEntry

from mp_api.client import MPRester

from rxn_network.enumerators.basic import BasicEnumerator
from rxn_network.enumerators.basic import BasicOpenEnumerator
from rxn_network.reactions.hull import InterfaceReactionHull
from collections.abc import Iterable

from pydmclab.data.thermochem import gas_thermo_data
from pydmclab.core.comp import CompTools

from pymatgen.entries.computed_entries import CompositionEnergyAdjustment
import math
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.entries.computed_entries import ComputedEntry
from solidstatesynth.gen.build_entry import BuildGibbsEntrySet,BuildGibbsEntry
from solidstatesynth.analyze.compound import AnalyzeChemsys
from rxn_network.costs.base import CostFunction

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
        stability_filter: float = 0.05
    ):
        """
        Initialize the calculator with the specified elements and temperature
        Args:
            precursors (Iterable[str]): Formula Strings of all the precursors that can be used in the reactions; precursors must span the whole chemical system besides Oxygen.
            targets (Iterable[str]): Formula Strings of all the targets to generate reactions for in the chemical system.
        """
        self.els = els
        if 'O' in self.els:
            self.els.extend(['H','C'])
        self._temperature = temperature
        self.with_theoretical = with_theoretical
        self.stability_filter = stability_filter
        self.precursors = AnalyzeChemsys(self.els).possible_precursors()
        # self.precursors = 
        
        #Find the chemical system to get query the MP entries from
        # elements = list(set([i for j in self._precursors+self._target+["O2"] for i in Composition(j).chemical_system_set]))
        # elements.sort()
        # self._chemsys = "-".join(elements)

        #Initialize the entries and reactions to None
        self.entries = None
        self.rxns = None
        self._build_calculator()


    def _get_entries(self) -> GibbsEntrySet:
        els = self.els
        temperature = self._temperature
        with_theoretical = self.with_theoretical
        stability_filter = self.stability_filter
        entry_set = BuildGibbsEntrySet(els = els,with_theoretical=with_theoretical,stability_filter=stability_filter).build_entry_set(temperature)
        return entry_set



    def _build_calculator(self):
        """
        Build the calculator by getting the entries and enumerating the reactions
        """
        self.entries = self._get_entries()
        print('entries obtained')
        #Use the BasicEnumerator and BasicOpenEnumerator to enumerate all the reactions from the precursors
        # self.rxns = BasicEnumerator(precursors=self._precursors, exclusive_precursors=True).enumerate(self.entries)
        self.rxns = BasicEnumerator().enumerate(self.entries)

        print('rxns enumerated')
        precursors = self.precursors
        self.rxns = self.rxns.add_rxn_set(BasicOpenEnumerator(open_phases=["O2"]).enumerate(self.entries))

        # self.rxns = self.rxns.add_rxn_set(BasicOpenEnumerator(open_phases=["O2"], precursors=self._precursors, exclusive_precursors=True).enumerate(self.entries))
        #Filter out duplicate reactions
        self.rxns = self.rxns.filter_duplicates()
        print('rxn duplicates filtered')
        return 


class TargetRxns():
    """
    class for calculating competition
    """
    def __init__(self, precursors, target, reactions, temperature = 300, environment = "air"):
        self.precursors = precursors
        self.target = target
        self.reactions = reactions
        self.temperature = temperature
        self.environment = environment
        self.NIST_gases = ['CO2','H2O']


    def _find_target_rxns(
            self   
        ) -> list[ComputedReaction]:
        """
        Helper method to find all the reactions that make the specified target without any byproducts (besides balancing gases)
        Args:
            target (str): String of the target to find the reactions for
            reactions (ReactionSet): ReactionSet of all the reactions to search through at the reaction temperature
        Returns:
            list[ComputedReaction]: List of all the reactions that make the target without any byproducts
        """
        reactions = self.reactions
        target = self.target
        
        target_rxns = []

        #The products that are allowed in the target reactions - the target and any balancing gases
        allowed_products = {Composition(target)} | set([Composition(i) for i in GAS_PARTIAL_PRESSURES["air"].keys()])

        #Check each reaction to see if its products are within the allowed products
        #The precursors should already be the only reactants in the reaction - besides O2 in some cases
        #These are already filtered by the correct atmosphere - if it is open to air, O2 can be a reactant

        for rxn in reactions:
            if set([i.reduced_composition for i in rxn.products]) <= allowed_products and (len(rxn.reactants) == 2 or len(set(rxn.reactants) - {Composition("O2")}) == 2):
                target_rxns.append(rxn)
        print('target rxns found')
        print(type(target_rxns))

        return ReactionSet.from_rxns(target_rxns)

    def get_rxns_with_desired_precursors(self):
        precursors = self.precursors
        environment = self.environment
        rxns = self._find_target_rxns()
        if precursors:
            if environment == "inert":
                new_rxns = rxns.get_rxns_by_reactants(reactants=precursors, return_set=True)
            else:
                self.precursors.append("O2")
                new_rxns = rxns.get_rxns_by_reactants(reactants=precursors, return_set=True)
        return new_rxns
    

    def _get_environment_correction(
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
        if environment not in GAS_PARTIAL_PRESSURES:
            return None

        #Composition for the gas
        comp_formula = CompTools(formula)
        # print(math.log(GAS_PARTIAL_PRESSURES[environment][formula]))
        #Calculate the environment correction for the gas on a per atom basis
        adjustment_per_atom = KB * math.log(GAS_PARTIAL_PRESSURES[environment][formula])*float(temperature)/comp_formula.n_atoms
        
        #Set a useful name for the correction - This is important because the correction does not automatically chenge with temperature
        name = f"{formula} {environment} correction @ {temperature}K"
        
        #Return the CompositionEnergyAdjustment object
        return CompositionEnergyAdjustment(adj_per_atom=adjustment_per_atom, n_atoms=comp_formula.n_atoms, name=name)


    def get_rxns_at_temp_env(
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
        rxns = self.get_rxns_with_desired_precursors()
            #Chemical potential of O at specified temperature and evironment partial pressure
        mu = KB * (temperature) * math.log(GAS_PARTIAL_PRESSURES[environment]["O2"])

        #Set the reactions to the new temperature
        #For inert environment, remove all reactions that contain O2 as a reactant
        #Determine the allowed gas products to add the correction to
        gas_comps = [i.composition for i in rxns.entries if i.reduced_formula in set(GAS_PARTIAL_PRESSURES["air"].keys())-{"O2"}]
        #Get the environment correction for each gas
        gas_dict = {i: self._get_environment_correction(i.reduced_formula, temperature = temp_new, environment = environment) for i in gas_comps}
        
        #Add the environment correction to the gas products
        for i in rxns.entries:
            if Composition(i.composition) in gas_comps:
                i.energy_adjustments.append(gas_dict[Composition(i.composition)])

        #Finally set O as an open element with the calculated chemical potential from partial pressure
        return rxns.set_chempot(open_el="O", chempot=mu)

    def metrics_at_temp_env(
        self,
        rxns = None
    ) -> list[dict[str, ComputedReaction | float]]:
        """
        Calculate the competition metrics for all the reactions at the specified temperature and environment
        Args:
            temp (float): Temperature in Kelvin to calculate the competition metrics at
            environment (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            metrics (list[dict[str, ComputedReaction | float]]): List of dictionaries with the reaction, energy, primary competition, secondary competition, and gamma cost function for each reaction
        """
        temperature = self.temperature
        environment = self.environment

        #Get the reactions at the correct temperature
        if not rxns:
            print('generating')
            rxns_at_temp = self.get_rxns_at_temp_env()
        else:
            print('using')
            rxns_at_temp = rxns
        print('rxns', rxns_at_temp)
        

        #Initialize the metrics data
        metrics = []

        #Loop through all the targets to calculate metrics for all the precursors combinations to make each target.
        # for target in self._targets:
            #Find the reactions that make the target without any byproducts - will calculate metrics for each reaction
        print(len(rxns_at_temp))
        for rxn in rxns_at_temp:
            print('rxn',rxn)
            #Get the precursors for the reaction - besides O2 if there are 2 solid precursors
            rxn_prec = set([i.reduced_composition for i in rxn.reactants])
            if len(rxn_prec) == 3 and Composition("O2") in rxn_prec:
                rxn_prec.remove(Composition("O2"))
            print(rxn_prec)
            #Get the filtered reactions for the InterfaceReactionHull that only contain the precursors in this reaction
            #In air, O2 is allowed as a reactant
            #In excess inert, O2 is not allowed as a reactant
            if environment == "air":
                filtered_rxns = list(rxns_at_temp.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]+["O2"]))
            else:
                filtered_rxns = list(rxns_at_temp.get_rxns_by_reactants([i.reduced_formula for i in rxn_prec]))

            #Build the InterfaceReactionHull from the precursors and filtered reactions
            irh = InterfaceReactionHull(*rxn_prec, filtered_rxns)

            #Calculate the competition metrics for the reaction and store them
            rxn_data = {}
            rxn_data["rxn"] = str(rxn)
            rxn_data["energy"] = rxn.energy_per_atom
            rxn_data["c1"] = irh.get_primary_competition(rxn)
            rxn_data["c2"] = irh.get_secondary_competition(rxn)
            rxn_data["gamma"] = 0.1 * rxn_data["energy"] + 0.45 * rxn_data["c1"] + 0.45 * rxn_data["c2"]
            
            #Useful for debugging:
            #rxn_data["competing_rxns"] = [str(i)+" dG="+str(i.energy_per_atom) for i in filtered_rxns]

            metrics.append(rxn_data)


        return metrics

    def get_rxns_at_new_temp_env(self, temp_new,rxns = None):
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
        temperature = self.temperature
        if not rxns:
            rxns = self.get_rxns_with_desired_precursors()
        else:
            rxns = rxns
        # precursors = self.precursors
        entries = rxns.entries
        indices = rxns.indices
        coeffs = rxns.coeffs
        gibbs_entries_at_temp = [i.get_new_temperature(temp_new) for i in entries if i.reduced_formula not in self.NIST_gases]
        nist_entries = [i for i in entries if i.reduced_formula in self.NIST_gases]
        nist_at_temp = [BuildGibbsEntry(i.reduced_formula, mp_data = None).gas_ExperimentalReferenceEntry_at_temp(temp_new) for i in nist_entries]
        gibbs_entries_at_temp.extend(nist_at_temp)
        new_rxns = ReactionSet(entries = gibbs_entries_at_temp, indices = indices, coeffs = coeffs)
        print(new_rxns)
        # new_rxns = ReactionSet.from_rxns(new_rxns)
            #Chemical potential of O at specified temperature and evironment partial pressure
        mu = KB * (temp_new) * math.log(GAS_PARTIAL_PRESSURES[environment]["O2"])

        #Set the reactions to the new temperature

        #For inert environment, remove all reactions that contain O2 as a reactant
        # if precursors:
        #     if environment == "inert":
        #         new_rxns = new_rxns.get_rxns_by_reactants(reactants=precursors, return_set=True)
        #     else:
        #         self.precursors.append("O2")
        #         new_rxns = new_rxns.get_rxns_by_reactants(reactants=precursors, return_set=True)
        # print('rxns by reactant', len(new_rxns))

        #Determine the allowed gas products to add the correction to
        gas_comps = [i.composition for i in new_rxns.entries if i.reduced_formula in set(GAS_PARTIAL_PRESSURES["air"].keys())-{"O2"}]
        
        #Get the environment correction for each gas
        gas_dict = {i: self._get_environment_correction(i.reduced_formula, temperature = temp_new, environment = environment) for i in gas_comps}
        
        #Add the environment correction to the gas products
        for i in new_rxns.entries:
            if Composition(i.composition) in gas_comps:
                i.energy_adjustments.append(gas_dict[Composition(i.composition)])

        #Finally set O as an open element with the calculated chemical potential from partial pressure
        return new_rxns.set_chempot(open_el="O", chempot=mu)


    # @property
    # def precursors(self) -> list[str]:
    #     """Precursors in the chemical system for the competition metrics"""
    #     return self._precursors
    
    # @property
    # def targets(self) -> list[str]:
    #     """Targets in the chemical system for the competition metrics"""
    #     return self._target
    
    # @property
    # def chemsys(self) -> str:
    #     """Chemical system for the competition metrics"""
    #     return self._chemsys
    
    

    

# def get_rxns_at_temp_env(
#         self,
#         temp_new = None,
#         no_byproducts = True,
#     ) -> ReactionSet:
#         """
