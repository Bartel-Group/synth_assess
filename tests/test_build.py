import unittest
from pydmclab.core.comp import CompTools
from pymatgen.core.composition import Composition
from rxn_network.entries.gibbs import GibbsComputedEntry
from solidstatesynth.gen.entries import Gibbs, GibbsSet
from pymatgen.ext.matproj import MPRester

# class UnitTestExtractData(unittest.TestCase):
#     sample_data = {'BaTiO3':}

class UnitTestBuildEntry(unittest.TestCase):

    def test_build_entry(self):
        formulas = ['BaTiO3','CO2','BaCO3']
        solids_data = {'BaTiO3':{'formation_energy_per_atom':  -3.4930003639166665,
                                 'volume':67.97032701923261, 'nsites':5},
                       'BaCO3':{'formation_energy_per_atom':  -2.6840060677500004, 
                                'volume':  319.7114081445892, 'nsites': 20}} 
        formula_energies = {'BaTiO3':{300:-3.322, 1450:-2.68},
                            'CO2':{300:-4.0876, 1450:-4.10700108825206},
                            'BaCO3':{300:-2.357, 1450:-1.789}}
        # temperatures = [300,1450]
        # environments = ['air','inert']
        # formula = formulas[0]
        for formula in formulas:
            n_atoms = CompTools(formula).n_atoms
            build_300 = Gibbs(formula, solids_data=solids_data, temperature=300)
            build_1450 = Gibbs(formula, solids_data=solids_data, temperature=1450)
            energy_300 = build_300.entry.energy/n_atoms
            energy_1450 = build_1450.entry.energy/n_atoms
            self.assertAlmostEqual(energy_300,formula_energies[formula][300], places=3)
            self.assertAlmostEqual(energy_1450,formula_energies[formula][1450], places=3)

    def test_carbonate_correction(self):
        formulas = ['CaCl2','CaCO3','CaMg(CO3)2']
