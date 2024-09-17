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

from pymatgen.entries.computed_entries import CompositionEnergyAdjustment
import math
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.entries.computed_entries import ComputedEntry

#Gas partial pressures in atm for different environments
PA_TO_ATM_CONV = 101300

GAS_PARTIAL_PRESSURES = {
    "air": {
        "O2": 21200 / PA_TO_ATM_CONV,
        "N2": 79033 / PA_TO_ATM_CONV,
        "CO2": 40.50 / PA_TO_ATM_CONV,
        #check with Nathan
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
#check with chris about inert values


#boltzmann constant
KB = 8.6173303e-5


class MetricsCalculator():
    """
    Class for calculating competition metrics for single-step reactions in a single chemical system.
    Uses the reaction-networks package for calculating reaction energy, primary and secondary competition, and the gamma cost function.
    Can be utilized by imnputting select precursors, targets, and temperature, and will output the competition metrics for each reaction.

    INSERT CREDIT TO MATT MCDERMOTT + maybe others???
    
    """

    def __init__(
        self,
        precursors: Iterable[str],
        targets: Iterable[str]
    ):
        """
        Initialize the calculator with the specified precursors and targets
        Args:
            precursors (Iterable[str]): Formula Strings of all the precursors that can be used in the reactions; precursors must span the whole chemical system besides Oxygen.
            targets (Iterable[str]): Formula Strings of all the targets to generate reactions for in the chemical system.
        """
        self._targets = [Composition(i).reduced_formula for i in targets]
        if precursors:
            print(precursors)
            self._precursors = [Composition(i).reduced_formula for i in precursors]
            elements = list(set([i for j in self._precursors+self._targets+["O2"] for i in Composition(j).chemical_system_set]))

        else:
            self._precursors = []
            elements = list(set([i for j in self._targets+["O2"] for i in Composition(j).chemical_system_set]))

        elements.sort()
        self._chemsys = "-".join(elements)

        #Initialize the entries and reactions to None
        self.entries = None
        self.rxns = None


        self._build_calculator()
        return




    def _build_calculator(self) -> None:
        """
        Queries MP for all the entries in the chemical system, and enumerates all the reactions from the precursors.
        All the reactions and entries are at 300K
        """

        #Query MP for all the ComputedStructureEntries in the chemical system
        with MPRester() as mpr:
            mp_entries = mpr.get_entries_in_chemsys(self._chemsys)


        #Create a GibbsEntrySet from the ComputedStructureEntries @ 300K
        self.entries = GibbsEntrySet.from_computed_entries(mp_entries, temperature=300, include_nist_data=True, apply_atmospheric_co2_correction=False)
        
        
        #Filter out any unstable phases 50 meV/atom above the hull
        self.entries = self.entries.filter_by_stability(0.05)        

        #Use the BasicEnumerator and BasicOpenEnumerator to enumerate all the reactions from the precursors
        self.rxns = BasicEnumerator(precursors=self._precursors).enumerate(self.entries)
        self.rxns = self.rxns.add_rxn_set(BasicOpenEnumerator(open_phases=["O2"], precursors=self._precursors).enumerate(self.entries))

        #Filter out duplicate reactions
        self.rxns = self.rxns.filter_duplicates()

        return

    def calculate_metrics_at_temp_env(
        self,
        temp: float = 1073,
        env: str = "air"
    ) -> list[dict[str, ComputedReaction | float]]:
        """
        Calculate the competition metrics for all the reactions at the specified temperature and environment
        Args:
            temp (float): Temperature in Kelvin to calculate the competition metrics at
            env (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            metrics (list[dict[str, ComputedReaction | float]]): List of dictionaries with the reaction, energy, primary competition, secondary competition, and gamma cost function for each reaction
        """

        #Get the reactions at the correct temperature
        rxns_at_temp = self.get_rxns_at_temp_env(temp=temp, env=env)

        #Initialize the metrics data
        metrics = []

        #Loop through all the targets to calculate metrics for all the precursors combinations to make each target.
        for target in self._targets:
            #Find the reactions that make the target without any byproducts - will calculate metrics for each reaction
            target_rxns = self._find_target_rxns(target, rxns_at_temp)

            for rxn in target_rxns:

                #Get the precursors for the reaction - besides O2 if there are 2 solid precursors
                rxn_prec = set([i.reduced_composition for i in rxn.reactants])
                if len(rxn_prec) == 3 and Composition("O2") in rxn_prec:
                    rxn_prec.remove(Composition("O2"))

                #Get the filtered reactions for the InterfaceReactionHull that only contain the precursors in this reaction
                #In air, O2 is allowed as a reactant
                #In excess inert, O2 is not allowed as a reactant
                
                if env == "air":
                    filtered_rxns = list(rxns_at_temp.get_rxns_by_reactants([i.formula for i in rxn_prec]+["O2"]))
                else:
                    filtered_rxns = list(rxns_at_temp.get_rxns_by_reactants([i.formula for i in rxn_prec]))

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





    def _find_target_rxns(
            self,
            target: str,
            reactions: ReactionSet,     
        ) -> list[ComputedReaction]:
        """
        Helper method to find all the reactions that make the specified target without any byproducts (besides balancing gases)
        Args:
            target (str): String of the target to find the reactions for
            reactions (ReactionSet): ReactionSet of all the reactions to search through at the reaction temperature
        Returns:
            list[ComputedReaction]: List of all the reactions that make the target without any byproducts
        """

        target_rxns = []

        #The products that are allowed in the target reactions - the target and any balancing gases
        allowed_products = {Composition(target)} | set([Composition(i) for i in GAS_PARTIAL_PRESSURES["air"].keys()])

        #Check each reaction to see if its products are within the allowed products
        #The precursors should already be the only reactants in the reaction - besides O2 in some cases
        #These are already filtered by the correct atmosphere - if it is open to air, O2 can be a reactant

        for rxn in reactions:
            if set([i.reduced_composition for i in rxn.products]) <= allowed_products:
                target_rxns.append(rxn)

        return target_rxns
    

    def _get_enironment_correction(
            self,
            formula: str,
            temp: float,
            environment: str = "air"
        ) -> CompositionEnergyAdjustment:
        """
        Get the environment correction for a gas at the specified temperature
        Different environments will have different partial pressures for gasses, which will affect the correction.
        Args:
            formula (str): Formula of the gas to get the environment correction for
            temp (float): Temperature in Kelvin to get the environment correction for
            environment (str): Environment to correct the entry for, either "air" for air environment or "inert" for excess inert environment
        Returns:
            CompositionEnergyAdjustment: Environment correction
        """
        if environment not in GAS_PARTIAL_PRESSURES:
            return None

        #Composition for the gas
        comp = Composition(formula)
    
        #Calculate the environment correction for the gas on a per atom basis
        adjustment_per_atom = KB * temp * math.log(GAS_PARTIAL_PRESSURES[environment][comp.reduced_formula]) / comp.num_atoms
        
        #Set a useful name for the correction - This is important because the correction does not automatically chenge with temperature
        name = f"{comp.reduced_formula} {environment} correction @ {temp}K"
        
        #Return the CompositionEnergyAdjustment object
        return CompositionEnergyAdjustment(adj_per_atom=adjustment_per_atom, n_atoms=comp.num_atoms, name=name)



    def get_entries_at_temp_env(
        self,
        temp: float = 800,
        env: str = "air"
    ) -> GibbsEntrySet:
        """
        Get the entries at the specified temperature and environment.
        This will set the correct partial pressure energy correction for NON-OXYGEN gases in the environment.
        NOTE: This does not correct the oxygen chemical potential, which is normally set during reaction enumeration.
        Args:
            temp (float): Temperature in Kelvin to get the entries at
            env (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            GibbsEntrySet: Entries at the specified temperature and environment
        """
        #Set the temperature for the new entries
        new_entries = self.entries.get_entries_with_new_temperature(new_temperature=temp)

        #Determine the allowed gas entries to add the correction to
        gas_comps = [i.composition for i in self.entries if i.reduced_formula in set(GAS_PARTIAL_PRESSURES["air"].keys())-{"O2"}]
        
        #Get the environment correction for each gas
        gas_dict = {i: self._get_enironment_correction(i.reduced_formula, temp, env) for i in gas_comps}

        #Apply the environment correction to each gas entry and save in new list
        corrected_gas_entries = [new_entries.get_adjusted_entry(new_entries.get_min_entry_by_formula(i.reduced_formula), gas_dict[i]) for i in gas_comps]

        #Remove the old gas entries the entry set
        for i in gas_comps:
            new_entries.discard(new_entries.get_min_entry_by_formula(i.reduced_formula))
        
        #Add the corrected gas entries to the entry set
        new_entries.update(corrected_gas_entries)

        return new_entries


    def get_rxns_at_temp_env(
        self,
        temp: float = 800,
        env: str = "air"
    ) -> ReactionSet:
        """
        Get the reactions at the specified temperature and environment. This requires adding a correction for the partial pressure of gases in the environment.
        For an air environment, O2 is a viable reactant, but not for an excess inert atmosphere.
        Args:
            temp (float): Temperature in Kelvin to get the reactions at
            env (str): Environment to correct the gas entries for ("air" or "inert")
        Returns:
            ReactionSet: Reactions at the specified temperature
        """
        #Chemical potential of O at specified temperature and evironment partial pressure
        mu = KB * (temp) * math.log(GAS_PARTIAL_PRESSURES[env]["O2"])

        #Set the reactions to the new temperature
        new_rxns = self.rxns.set_new_temperature(new_temp=temp)

        #For inert environment, remove all reactions that contain O2 as a reactant
        if env == "inert":
            new_rxns = new_rxns.get_rxns_by_reactants(reactants=self._precursors, return_set=True)
        

        #Determine the allowed gas products to add the correction to
        gas_comps = [i.composition for i in self.entries if i.reduced_formula in set(GAS_PARTIAL_PRESSURES["air"].keys())-{"O2"}]
        
        #Get the environment correction for each gas
        gas_dict = {i: self._get_enironment_correction(i.reduced_formula, temp, env) for i in gas_comps}
        
        #Add the environment correction to the gas products
        for i in new_rxns.entries:
            if Composition(i.composition) in gas_comps:
                i.energy_adjustments.append(gas_dict[Composition(i.composition)])

        #Finally set O as an open element with the calculated chemical potential from partial pressure
        return new_rxns.set_chempot(open_el="O", chempot=mu)


    @property
    def precursors(self) -> list[str]:
        """Precursors in the chemical system for the competition metrics"""
        return self._precursors
    
    @property
    def targets(self) -> list[str]:
        """Targets in the chemical system for the competition metrics"""
        return self._targets
    
    @property
    def chemsys(self) -> str:
        """Chemical system for the competition metrics"""
        return self._chemsys
    
    

    