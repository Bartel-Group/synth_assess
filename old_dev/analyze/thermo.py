
import os
# from solidstatesynth.data.tm import (
#     get_tm_precursors,
#     get_tm_targets,
#     get_mp_cmpds,
# )
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json
from itertools import combinations
from pymatgen.ext.matproj import MPRester
from pymatgen.entries.computed_entries import GibbsComputedStructureEntry
from rxn_network.entries.gibbs import GibbsComputedEntry
from rxn_network.entries.nist import NISTReferenceEntry
from pymatgen.core.composition import Composition

DATADIR = "../data"

class AnalyzeThermo(object):
    def init(self):
        self.tm_precursors = read_json(os.path.join(DATADIR, 'tm_precursors.json'))['data']
        self.tm_targets = read_json(os.path.join(DATADIR, 'tm_targets.json'))['data']
        self.mp_cmpds = read_json(os.path.join(DATADIR, 'mp_cmpds.json'))['data']
        self.mp_icsd_cmpds = read_json(os.path.join(DATADIR, 'mp_icsd_cmpds.json'))['data']

    # @property
    def mp_data_from_formula(self,formula):
        """
        Returns:
            MP energy and volume associated with ground state of the compound
        """
        data = MPRester().materials.summary.search(formula = formula, fields=["energy_per_atom","volume"])
        return data
    
    # @property
    def ground_GibbsComputedEntry_at_temp(self,formula, temp, stability_filter = None):
        """
        Args:
            formula (str) : formula
            stability_filter (float) : energy per atom threshold for stability (or None if filter not desired)
        Returns:
            GibbsComputedStructureEntry for ground state polymorph at 300 K at a temperature of interest
        """
        mp_data = self.mp_data_from_formula(formula)
        # gibbs_computed_structure_entries = [GibbsComputedStructureEntry(
        #                                 structure = a.structure, 
        #                                 formation_enthalpy_per_atom = a.energy_per_atom, 
        #                                 composition=Composition(formula),
        #                                 temp = 300) 
        #                                 for a in mp_data]
        # gd_state_entry = gibbs_computed_structure_entries[0]
        # for entry in gibbs_computed_structure_entries:
        gibbs_computed_entries = [GibbsComputedEntry(
                                        volume_per_atom = (a.volume)/(CompTools(formula).n_atoms), 
                                        formation_energy_per_atom = a.energy_per_atom, 
                                        composition=Composition(formula),
                                        temperature = 300) 
                                        for a in mp_data]
        gd_state_entry = gibbs_computed_entries[0]
        for entry in gibbs_computed_entries:
            if entry.energy_per_atom < gd_state_entry.energy_per_atom:
                gd_state_entry = entry
        if stability_filter:
            if gd_state_entry.energy_per_atom > stability_filter:
                return None
        if temp == 300:
            return gd_state_entry
        print('volume per atom', formula,gd_state_entry.volume_per_atom)
        return gd_state_entry.get_new_temperature(temp)
    
    def gas_GibbsComputedEntry_at_temp(self,formula, temp = 300):
        return NISTReferenceEntry(composition=Composition(formula), temperature =temp)
    
    
        

def main():
    return


if __name__ == "__main__":
    main()