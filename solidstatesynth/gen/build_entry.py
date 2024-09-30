
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
            True if the target is a NIST gas
        """
        return True if self.formula in ['CO2','H2O','C1O2','H2O1'] else False

    # def get_temperature_interpolated_entry(self,temp,energies):
    #     upper_T = math.ceil(temp / 100) * 100
    #     lower_T = upper_T - 100
    #     if temp >= 2000:
    #         lower_T = 1900
    #         upper_T = 2000
    #     if temp <= 300:
    #         lower_T = 300
    #         upper_T = 400

    #     temperatures_to_grab = [str(lower_T), str(upper_T)]
    #     dGfs = [energies[t]['Ef'] for t in temperatures_to_grab]
    #     if not dGfs:
    #         return None
    #     if dGfs[0] == None or dGfs[1] == None:
    #         return None
    #     slope = (dGfs[1] - dGfs[0]) / (
    #         float(temperatures_to_grab[1]) - float(temperatures_to_grab[0])
    #     )
    #     # interpolation-- difference in dGf divided by difference in temperature
    #     return dGfs[0] + slope * (temp - float(temperatures_to_grab[0]))
        
    def gas_ExperimentalReferenceEntry_at_temp(self, temp = 300):
        formula = self.formula
        refs = self.gas_data
        ExperimentalReferenceEntry.REFERENCES = refs
        # print('refs',ExperimentalReferenceEntry.REFERENCES)
        exp_entry = ExperimentalReferenceEntry(composition=Composition(formula), temperature = temp)
        return exp_entry

    # def gas_ComputedEntry_at_temp(self, temp):
    #     """
    #     Returns:
    #         GibbsComputedStructureEntry for H2O at 300 K
    #     """
    #     formula = self.formula
    #     energies = self.gas_data[self.formula]
    #     dGf = self.get_temperature_interpolated_entry(temp,energies)
    #     return ComputedEntry(composition=Composition(formula), energy = dGf)
    
    # @property
    def ground_GibbsComputedEntry_at_temp(self,formula, temp=300):
        """
        Args:
            formula (str) : formula
        Returns:
            GibbsComputedStructureEntry for ground state polymorph at 300 K at a temperature of interest
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
    def __init__(self,target, with_theoretical = False, stability_filter = 0.05):
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
        self.target_els = target_els

    def is_competing_formula(self,mp_formula):
        
        """
        mp_entry (dict) : entry from MP data of the form {'formula': str, 'energy_above_hull': float,...}
        """
        target_els = self.target_els
        formula_els = CompTools(mp_formula).els
        if all([el in target_els for el in formula_els]):
            return True
        return False

    def target_competing_formulas(self):
    
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