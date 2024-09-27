
from pydmclab.utils.handy import read_json,write_json 
import os
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.entry_set import GibbsEntrySet
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.core.composition import Composition
from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeTarget
from pydmclab.core.comp import CompTools
from solidstatesynth.extract.mp import get_mp_experimental

DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"


class BuildGibbsEntry:
    def __init__(self,formula, mp_data):
        self.formula = CompTools(formula).clean
        self.mp_data = mp_data

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
            True if the target is a NIST gas
        """
        return True if self.formula in ['C1O2','H2O1'] else False
    
    # @property
    def ground_GibbsComputedEntry_at_temp(self,formula, temp=300):
        """
        Args:
            formula (str) : formula
            stability_filter (float) : energy per atom threshold for stability (or None if filter not desired)
        Returns:
            GibbsComputedStructureEntry for ground state polymorph at 300 K at a temperature of interest
        """
        if type(formula) == dict:
            mp_data = formula
        else:
            mp_data = self.mp_data_for_formula(formula)
        volume_per_atom_calc = mp_data['volume']/mp_data['nsites']
        gibbs_computed_entry = GibbsComputedEntry(
                                        volume_per_atom = volume_per_atom_calc, 
                                        formation_energy_per_atom = mp_data['formation_energy_per_atom'], 
                                        composition=Composition(formula),
                                        temperature = temp) 
        return gibbs_computed_entry
    
    def gas_GibbsComputedEntry_at_temp(self,formula, temp = 300):
        if type(formula) == dict:
            formula = formula['formula']
        return NISTReferenceEntry(composition=Composition(formula), temperature = temp)
    
class BuildGibbsEntrySet():
    def __init__(self,target, with_theoretical = False, stability_filter = 0.1):
        self.target = target
        self.with_theoretical = with_theoretical
        self.stability_filter = stability_filter
        if not with_theoretical:
            self.mp_data = read_json(os.path.join(DATADIR, '240926_mp_experimental.json'))['data']
        else:
            self.mp_data = read_json(os.path.join(DATADIR, '240925_mp_ground_data.json'))['data']

    def is_competing_formula(self,mp_entry):
        
        """
        mp_entry (dict) : entry from MP data of the form {'formula': str, 'energy_above_hull': float,...}
        """
        target = self.target
        stability_filter = self.stability_filter
        target_els = CompTools(target).els
        target_els.extend(AnalyzeTarget(target).flexible_els)
        formula_els = CompTools(mp_entry['formula']).els
        formula_stability = mp_entry['energy_above_hull']
        if all([el in target_els for el in formula_els]):
            if formula_stability < stability_filter:
                return True
        return False
    
    def target_competing_formulas(self):
        mp_data = self.mp_data
        mp_comp_entries = [mp_data[formula] for formula in mp_data if self.is_competing_formula(mp_data[formula])]
        return mp_comp_entries

    def build_entry_set(self, temp):
        competing_formulas = self.target_competing_formulas()
        mp_data = self.mp_data
        entries = []
        for formula in competing_formulas:
            GibbsEntry = BuildGibbsEntry(formula, mp_data)
            if GibbsEntry.is_NIST_gas:
                entry = GibbsEntry.gas_GibbsComputedEntry_at_temp(formula, temp)
            else:
                entry = GibbsEntry.ground_GibbsComputedEntry_at_temp(formula, temp)
            entries.append(entry)
        return GibbsEntrySet(entries)

def main():
    return

if __name__ == '__main__':
    main()