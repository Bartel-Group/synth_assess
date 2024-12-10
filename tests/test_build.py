import unittest
from pydmclab.core.comp import CompTools
from pymatgen.core.composition import Composition
from rxn_network.entries.gibbs import GibbsComputedEntry
from solidstatesynth.gen.entries import Gibbs, GibbsSet, FormulaChecker
from solidstatesynth.gen.metrics_calculation import EnumerateRxns
from pymatgen.ext.matproj import MPRester

# class UnitTestExtractData(unittest.TestCase):
#     sample_data = {'BaTiO3':}

class UnitTestGibbs(unittest.TestCase):

    def test_is_carbonate(self):
        formulas = ['BaTiO3','CO2','BaCO3']
        carbonates = {'BaTiO3':False, 'CO2':False, 'BaCO3':True}
        for formula in formulas:
            self.assertEqual(Gibbs(formula).is_carbonate, carbonates[formula])

    def test_carbonate_correction(self):
        formulas = ['CaCl2','CaCO3','CaMg(CO3)2']
        carbonate_corrections = {'CaCl2':0, 'CaCO3':0.830, 'CaMg(CO3)2':1.660}
        for formula in formulas:
            correction = Gibbs(formula).carbonate_correction
            self.assertEqual(correction,carbonate_corrections[formula])

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

class UnitTestFormulaChecker(unittest.TestCase):
    def test_sub_chemsys(self):
        formulas = ['BaTiO3','CO2','BaCO3']
        chemsys_els = ['Ba','Ti','O']
        subsys = {'BaTiO3':True, 'CO2':False, 'BaCO3':False}
        for formula in formulas:
            self.assertEqual(
                FormulaChecker(formula, chemsys_els).sub_chemsys, 
                subsys[formula]
                )


    def test_is_relevant(self):
        formulas = ['BaTiO3','CO2','BaCO3', 'CaCO3']
        chemsys_els = ['Ba','Ti','O']
        relevant = {'BaTiO3':True, 'CO2':True, 'BaCO3':True, 'CaCO3':False}
        for formula in formulas:
            self.assertEqual(
                FormulaChecker(formula, chemsys_els).is_relevant, 
                relevant[formula]
                )

class UnitTestGibbsSet(unittest.TestCase):
    def test_formulas(self):
        chemsys_els = [['Cs','Cl'],['Cs','O']]
        els_data = {'CsO2': {},'Cs4O3': {},'Cs7O': {},
                    'Cs3O': {},'CsO3': {},'Cs2O3': {},
                    'Cs2O2': {},'O2': {},'Cs2O': {},
                    'Cs11O3': {},'Cs': {},'CsOH': {},
                    'Cs2CO3': {},'CsCl': {},'Cl2': {},
                    'CO2': {},'H2O': {}}
        f_true = list(els_data.keys())
        entries = {'Cs-Cl': ['CsCl', 'Cl2', 'Cs'], 
                   'Cs-O': [e for e in f_true if 'Cl' not in e]}
        for entry in chemsys_els:
            chemsys = '-'.join(entry)
            gibbs_set = GibbsSet(entry, solids_data = els_data)
            self.assertEqual(list(set(gibbs_set.formulas)), list(set(entries[chemsys])))

class UnitTestEnumerate(unittest.TestCase):
    def test_get_entries(self):




