
from pydmclab.utils.handy import read_json,write_json 
import os
import math
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.core.composition import Composition
from pymatgen.entries.computed_entries import ComputedEntry
from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeTarget
from pydmclab.core.comp import CompTools
from pydmclab.data.thermochem import gas_thermo_data
from solidstatesynth.extract.mp import get_useful_mp_data, get_gases_data
from rxn_network.entries.experimental import ExperimentalReferenceEntry


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
        return mp_data[formula]
       
    
    @property
    def is_NIST_gas(self):
        """
        Returns:
            True if the target is a NIST gas (accounting for )
        """
        return True if self.formula in ['C1O2','H2O1'] else False
        
    def gas_ExperimentalReferenceEntry_at_temp(self, temp = 300):
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
        exp_entry = ExperimentalReferenceEntry(composition=Composition(formula), temperature = temp)
        return exp_entry

    def ground_GibbsComputedEntry_at_temp(self,formula, temp=300):
        """
        Args:
            formula (str) : formula
        Returns:
            GibbsComputedEntry for ground state polymorph at 300 K at a temperature of interest
        """
        if type(formula) == dict:
            mp_data = formula
        else:
            mp_data = self.mp_data_for_formula
        volume_per_atom_calc = mp_data['volume']/mp_data['nsites']
        gibbs_computed_entry = GibbsComputedEntry(
                                        volume_per_atom = volume_per_atom_calc, 
                                        formation_energy_per_atom = mp_data['formation_energy_per_atom'], 
                                        composition=Composition(formula),
                                        temperature = temp) 
        return gibbs_computed_entry
    
    
class BuildGibbsEntrySet():
    """
    Builds an EntrySet for formulas associated with a target of interest. This entry set is comprised
    of GibbsComputedEntries for ground state polymorphs and ExperimentalReferenceEntries for gases.
    Note that initialization requires a determination of whether hypothetical compounds may be accounted for 
    (with theoretical) and a stability filter (for energy above the hull)
    Depending on the with_theoretical initialization, MP with or without theoretical compounds may be employed 
    """
    def __init__(self,target, with_theoretical = True, stability_filter = 0.05):

        self.target = target
        self.with_theoretical = with_theoretical
        self.stability_filter = stability_filter
        if not with_theoretical:
            self.mp_data = read_json(os.path.join(DATADIR, '240926_mp_experimental.json'))['data']
        else:
            self.mp_data = read_json(os.path.join(DATADIR, '240925_mp_ground_data.json'))['data']
        mp_data_new = get_useful_mp_data(self.mp_data,n_els_max = CompTools(target).n_els)
        self.mp_formulas = list(mp_data_new.keys())
        target_els = CompTools(target).els
        target_els.extend(AnalyzeTarget(target).flexible_els)
        # should H, C be included for competing reactions? I would imagine yes?
        self.target_els = target_els

    def is_competing_formula(self,mp_formula):
        
        """
        Returns:
            True if the formula is in the chemical system (or sub chemical system) of the target
        """
        target_els = self.target_els
        formula_els = CompTools(mp_formula).els
        if all([el in target_els for el in formula_els]):
            return True
        return False

    def target_competing_formulas(self):
        """
        Returns a list of all relevant competing formulas from Materials Project for the target (as defined above)
        Note that depending on which version of MP is used, this may account only for ground state experimental formulas or all ground state formulas.
        """
    
        mp_formulas = self.mp_formulas
        formulas = []
        for formula in mp_formulas: 
            # print('formula',formula)
            if self.is_competing_formula(formula):
                print('competing formula')
                formulas.append(formula)
        return formulas

    def build_entry_set(self, temp):
        competing_formulas = self.target_competing_formulas()
        print('competing formulas',len(competing_formulas))
        mp_data = self.mp_data
        entries = []
        for formula in competing_formulas:
            print('formula',formula)
            GibbsEntry = BuildGibbsEntry(formula, mp_data)
            if GibbsEntry.is_NIST_gas:
                entry = GibbsEntry.gas_ExperimentalReferenceEntry_at_temp(temp)
                print(entry)
            else:
                entry = GibbsEntry.ground_GibbsComputedEntry_at_temp(formula, temp)
            entries.append(entry)
            print('entry added')
        return GibbsEntrySet(entries)

def main():
    return

if __name__ == '__main__':
    main()