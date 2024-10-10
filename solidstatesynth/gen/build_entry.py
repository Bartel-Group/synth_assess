
from pydmclab.utils.handy import read_json,write_json 
import os
import math
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from pydmclab.core.comp import CompTools
from pydmclab.data.thermochem import gas_thermo_data
from solidstatesynth.extract.mp import get_gases_data
from rxn_network.entries.experimental import ExperimentalReferenceEntry
from itertools import combinations


DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"


class BuildGibbsEntry:
    def __init__(self,formula, mp_data):
        self.formula = CompTools(formula).clean
        self.mp_data = mp_data
        self.gas_data = get_gases_data()

    @property
    def mp_data_for_formula(self):
        """
        Returns:
            MP energy and volume associated with ground state of the compound
        """
        mp_data = self.mp_data
        formula = self.formula
        formula_data = [entry for entry in mp_data if CompTools(entry['formula']).clean == formula]
        if formula_data:
            return formula_data[0]
        return None
       
    
    @property
    def is_NIST_gas(self):
        """
        Returns:
            True if the target is a NIST gas (accounting for )
        """
        return True if self.formula in ['C1O2','H2O1'] else False
        
    def gas_ExperimentalReferenceEntry_at_temp(self, temperature = 300):
        """
        Returns entry set for a gas at a temperature of interest. Note that the structure
        of this entry does not allow for modification at different temperatures (because
        temperatures are experimentally determined and G(T) cannot be extrapolated). The data
        used here in references is NIST experimental data for H2O and CO2 at different temperatures.
        NOTE THAT the functionality of this code requires that H2O and CO2 be written as such rather
        than clean formulas H2O1 and C1O2. The get_gases_data function in extract/mp.py is written
        such that keys are in this format as well. Note that temperature keys are ints (not strings)
        """
        formula = self.formula
        refs = self.gas_data
        ExperimentalReferenceEntry.REFERENCES = refs
        exp_entry = ExperimentalReferenceEntry(composition=Composition(formula), temperature = temperature)
        return exp_entry


    def ground_GibbsComputedEntry_at_temp(self,formula_data, temperature=300):
        """
        Args:
            formula (str) or formula MP entry : formula
        Returns:
            GibbsComputedEntry for ground state polymorph at 300 K at a temperature of interest
        """
        # if type(formula) == dict:
        # else:
        #     mp_data = self.mp_data_for_formula
        volume_per_atom_calc = formula_data['volume']/formula_data['nsites']
        gibbs_computed_entry = GibbsComputedEntry(
                                        volume_per_atom = volume_per_atom_calc, 
                                        formation_energy_per_atom = formula_data['formation_energy_per_atom'], 
                                        composition=Composition(formula_data['formula']),
                                        temperature = temperature) 
        return gibbs_computed_entry
    
    
class BuildGibbsEntrySet():
    """
    Builds an EntrySet for formulas associated with a target of interest. This entry set is comprised
    of GibbsComputedEntries for ground state polymorphs and ExperimentalReferenceEntries for gases.
    Note that initialization requires a determination of whether hypothetical compounds may be accounted for 
    (with theoretical) and a stability filter (for energy above the hull)Depending on the with_theoretical initialization, 
    MP with or without theoretical compounds may be employed 

    """
    def __init__(self,els, with_theoretical = True, stability_filter = 0.05, formulas = None):

        self.with_theoretical = with_theoretical
        self.stability_filter = stability_filter
        if formulas:
            self.formulas = formulas
        else:
            if not with_theoretical:
                self.data = read_json(os.path.join(DATADIR, '241002_mp_experimental.json'))['data']
            else:
                self.data = read_json(os.path.join(DATADIR, '241002_mp_gd.json'))['data']
            self.formulas = [entry['formula'] for entry in self.data]
        self.target_els = els

    def is_competing_formula(self,mp_formula):
        
        """
        Returns:
            True if the formula is in the chemical system (or sub chemical system) of the target
        """
        target_els = self.target_els
        allowed_els = []
        new_allowed_els = []
        for n in range(1, len(target_els)):
            allowed_els.extend(list(combinations(target_els, n)))
        if 'O' in target_els:
            flexible_els = ['C','H']
        for el in flexible_els:
            for el_combo in allowed_els:
                el_combo = list(el_combo)
                print(el_combo)
                if ("O" in el_combo) and (el not in el_combo) and (len(el_combo) > 1):
                    el_combo.append(el)
                    el_combo = tuple(sorted(el_combo))
                    new_allowed_els.append(el_combo)
        print(new_allowed_els)
        allowed_els = set(allowed_els + new_allowed_els)
        print(allowed_els)
        # filter our big list of precursors down to those that we deemed "possible"
        # precursors = [p for p in precursors if tuple(CompTools(p).els) in allowed_els]
        formula_els = CompTools(mp_formula).els
        if tuple(CompTools(mp_formula).els) in allowed_els:
                return True
        elif all([el in target_els for el in formula_els]):
            return True
        return False

    def chemsys_competing_formulas(self):
        """
        Returns a list of all relevant competing formulas from Materials Project for the target (as defined above)
        Note that depending on which version of MP is used, this may account only for ground state experimental formulas or all ground state formulas.
        """
        mp_data = self.data
        formulas_new = []
        formula_strings = []
        for entry in mp_data: 
            formula_str = entry['formula']
            if formula_str not in formula_strings:
                if self.is_competing_formula(formula_str):
                    formulas_new.append(entry)
                    formula_strings.append(formula_str)
        return formulas_new

    def build_entry_set(self, temperature=300):
        """
        Builds an entry set for the chemical space of interest. This entry set is comprised of GibbsComputedEntries for 
        ground state polymorphs and ExperimentalReferenceEntries for gases. Temperature can be specified.
        """
        competing_formulas = self.chemsys_competing_formulas()
        print('competing formulas found')
        data = self.data
        entries = []
        for formula_data in competing_formulas:
            formula_str = formula_data['formula']
            print('building entry')
            GibbsEntry = BuildGibbsEntry(formula_str, data)
            if GibbsEntry.is_NIST_gas:
                entry = GibbsEntry.gas_ExperimentalReferenceEntry_at_temp(temperature)
            else:
                entry = GibbsEntry.ground_GibbsComputedEntry_at_temp(formula_data, temperature)
            entries.append(entry)
            print('entry added')
        return GibbsEntrySet(entries)

def main():
    return

if __name__ == '__main__':
    main()