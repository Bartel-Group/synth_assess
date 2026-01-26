
import unittest
from pydmclab.core.comp import CompTools
from synth_assess.selectivity.entries import Gibbs, GibbsSet, FormulaChecker
from synth_assess.selectivity.rxn_metrics import PrecursorSet, EnumerateRxns, TempEnvCorrections, RxnsAtNewTempEnv, AnalyzeReactionSet


solids_data = {'Cs1O3': {'volume': 353.16244441729305,
  'formation_energy_per_atom': -0.6521186182327591,
  'nsites': 16,'tm_precursor': False,'theoretical': False},
 'Cs2O2': {'volume': 110.46533766915991,
  'formation_energy_per_atom': -1.230510951465517,
  'nsites': 4,'tm_precursor': False,'theoretical': False},
 'Cs2O1': {'volume': 115.60293111947884,
  'formation_energy_per_atom': -1.188630010287356,
  'nsites': 3,'tm_precursor': True,'theoretical': False},
 'Cs1': {'volume': 3393.9323239909295,
  'formation_energy_per_atom': 0.0,
  'nsites': 29,'tm_precursor': True,'theoretical': False},
 'C1Cs2O3': {'volume': 534.150691010004,
  'formation_energy_per_atom': -2.067810993477012,
  'nsites': 24,'tm_precursor': True,'theoretical': False},
 'Cl1': {'volume': 140.7080127118759,
  'formation_energy_per_atom': 0.0,
  'nsites': 4,'tm_precursor': True,'theoretical': False},
 'H2O1': {'volume': 124.89445215830659,
  'formation_energy_per_atom': -1.281629922499999,
  'nsites': 12,'tm_precursor': True,'theoretical': False},
 'C1Cd1O3': {'volume': 120.68862467740928,
  'formation_energy_per_atom': -1.762851043499999,
  'nsites': 10,'tm_precursor': True,'theoretical': False},
 'C1O3': {'volume': 32.99685994750084,
  'formation_energy_per_atom': 1.746732446874999,
  'nsites': 4,'tm_precursor': False,'theoretical': True},
 'Cd1O2': {'volume': 157.87700371639545,
  'formation_energy_per_atom': -0.8535735183333331,
  'nsites': 12,'tm_precursor': False,'theoretical': False},
 'Cd1': {'volume': 23.25485140168592,
  'formation_energy_per_atom': 0.0,
  'nsites': 1,'tm_precursor': True,'theoretical': False},
 'S1': {'volume': 1152.997853241063,
  'formation_energy_per_atom': 0.0,
  'nsites': 32,'tm_precursor': True,'theoretical': False},
 'Cd1S1': {'volume': 104.86287534018373,
  'formation_energy_per_atom': -0.8838398540624991,
  'nsites': 4,'tm_precursor': True,'theoretical': False},
 'Cd1S2': {'volume': 263.1219868050095,
  'formation_energy_per_atom': -0.659416383749999,
  'nsites': 12,'tm_precursor': False,'theoretical': False},
 'Cd2S1': {'volume': 78.99910856081462,
  'formation_energy_per_atom': -0.36680444437500004,
  'nsites': 3,'tm_precursor': False,'theoretical': True},
 'Cs2O3': {'volume': 483.0414253319752,
  'formation_energy_per_atom': -1.028256548672414,
  'nsites': 20,'tm_precursor': False,'theoretical': True},
 'O1': {'volume': 107.6438082346079,
  'formation_energy_per_atom': 0.0,'nsites': 8,
  'tm_precursor': True,'theoretical': False},
 'Cs1H1O1': {'volume': 123.32082021448794,
  'formation_energy_per_atom': -1.492687159310345,
  'nsites': 6,'tm_precursor': False,'theoretical': True},
 'Cl1Cs1': {'volume': 88.04926423740011,
  'formation_energy_per_atom': -2.254950123965517,
  'nsites': 2,'tm_precursor': True,'theoretical': False},
 'C1O2': {'volume': 195.38450985171852,
  'formation_energy_per_atom': -1.769495735833333,
  'nsites': 12,'tm_precursor': True,'theoretical': False}}
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

# helper function
def get_reaction_dict_from_string(reaction_string):
    """
    Args: reaction string
    Returns: dictionary with keys 'reactants' and 'products' where
    the values are lists of the reactants and products in the reaction
    *** IF CLEANABLE -- otherwise returns None ***
    Uses: this is the most usable reaction format from which to calculate
    dG_rxn using the pydmclab ReactionEnergy class-- get_dGrxn_at_T takes
    a reaction dictionary as an argument
    """
    reactant_list = []
    reactant_coeffs = []
    product_list = []
    product_coeffs = []
    if '->' in reaction_string:
        reactants, products = reaction_string.split(" -> ")
        # print(reactants)
    elif '==' in reaction_string:
        reactants, products = reaction_string.split(" == ")
    if "+" in reactants:
        reactants = reactants.split(" + ")
    else:
        reactants = [reactants]
    # print('reactants', reactants)
    for reactant in reactants:
        if " " in reactant:
            if len(reactant.split(" ")) != 1:
                coefficient, reactant = reactant.split(" ")
                if float(coefficient) > 0:
                    try:
                        CompTools(reactant).clean
                        reactant_list.append(reactant)
                        reactant_coeffs.append(coefficient)
                    except:
                        return None
        else:
            reactant_list.append(reactant)
            reactant_coeffs.append(1)
    if "+" in products:
        products = products.split(" + ")
    else:
        products = [products]
    for product in products:
        if " " in product:
            if len(product.split(" ")) != 1:
                coefficient, product = product.split(" ")
                if float(coefficient) > 0:
                    try:
                        CompTools(product).clean
                        product_list.append(product)
                        product_coeffs.append(coefficient)
                    except:
                        return None
        else:
            product_list.append(product)
            product_coeffs.append(1)
    # print(reactant_list,product_list)
    return {"reactants": reactant_list, "products": product_list, 
            "reactant_coeffs": reactant_coeffs, "product_coeffs": product_coeffs}

class UnitTestGibbs(unittest.TestCase):

    def test_is_carbonate(self):
        formulas = ['BaTiO3','CO2','BaCO3']
        carbonates = {'BaTiO3':False, 'CO2':False, 'BaCO3':True}
        for formula in formulas:
            self.assertEqual(Gibbs(formula, carbonates).is_carbonate, carbonates[formula])

    def test_carbonate_correction(self):
        formulas = ['CaCl2','CaCO3','CaMg(CO3)2']
        carbonate_corrections = {'CaCl2':0, 'CaCO3':0.166, 'CaMg(CO3)2':0.166}
        for formula in formulas:
            correction = Gibbs(formula, carbonate_corrections).carbonate_correction
            self.assertAlmostEqual(correction,carbonate_corrections[formula], delta=0.003)

    def test_build_entry(self):
        formulas = ['Cd1O2','C1O2','C1Cd1O3']
        formula_energies = {'Cd1O2':{300:-0.657, 1450:0.207},
                            'C1O2':{300:-4.0876, 1450:-4.10700108825206},
                            'C1Cd1O3':{300:-1.428, 1450:-0.757}}

        for formula in formulas:
            formula_new = CompTools(formula).clean
            n_atoms = CompTools(formula_new).n_atoms
            build_300 = Gibbs(formula_new, solids_data=solids_data, temperature=300)
            build_1450 = Gibbs(formula_new, solids_data=solids_data, temperature=1450)
            energy_300 = build_300.entry.energy
            energy_1450 = build_1450.entry.energy
            if CompTools(formula).clean != 'C1O2':
                energy_300 = energy_300/n_atoms
                energy_1450 = energy_1450/n_atoms
            self.assertAlmostEqual(energy_300,formula_energies[formula][300], delta= 0.003)
            self.assertAlmostEqual(energy_1450,formula_energies[formula][1450], delta = 0.003)

class UnitTestFormulaChecker(unittest.TestCase):
    def test_sub_chemsys(self):
        formulas = ['Ba1Ti1O3','C1O2','Ba1C1O3']
        chemsys_els = ['Ba','Ti','O']
        subsys = {'Ba1Ti1O3':True, 'C1O2':False, 'Ba1C1O3':False}
        for formula in formulas:
            self.assertEqual(
                FormulaChecker(formula, chemsys_els).sub_chemsys, 
                subsys[formula]
                )
            

    def test_is_relevant(self):
        formulas = ['Ba1Ti1O3','Ba1C1O3', 'Ca1C1O3']
        chemsys_els = ['Ba','Ti','O']
        relevant = {'Ba1Ti1O3':True, 'Ba1C1O3':True, 'Ca1C1O3':False}
        for formula in formulas:
            self.assertEqual(
                FormulaChecker(formula, chemsys_els).is_relevant, 
                relevant[formula]
                )

class UnitTestGibbsSet(unittest.TestCase):
    def test_formulas(self):
        chemsys_els = [['Cs','Cl'],['Cs','O']]
        entries = {'Cs-Cl': ['Cl1Cs1', 'Cl1', 'Cs1'], 
                   'Cs-O': ['Cs1O3','Cs2O3','Cs2O2','O1','Cs2O1',
                            'Cs1','Cs1H1O1','C1Cs2O3','C1O3', 'C1O2','H2O1']}
        for entry in chemsys_els:
            chemsys = '-'.join(entry)
            gibbs_set = GibbsSet(entry, solids_data = solids_data)
            gibbs_set_formulas = [CompTools(f).clean for f in gibbs_set.formulas]
            self.assertTrue(all(i in gibbs_set_formulas for i in entries[chemsys]))

class UnitTestPrecursorSet(unittest.TestCase):
    def test_precursors(self):
        els = ['Cs', 'O']
        prec_tm = set([CompTools(p).clean for p in PrecursorSet(els,solids_data).precursors])
        prec = set([CompTools(p).clean for p in PrecursorSet(els,solids_data, restrict_to_tm_precursors=False).precursors])
        print('prec', prec)
        true_prec_tm = set(['O1', 'Cs2O1', 'Cs1', 'C1Cs2O3'])
        true_prec = set(['Cs2O1','Cs1O3','Cs2O2', 'O1', 'Cs1', 'C1Cs2O3'])
        self.assertEqual(prec_tm, true_prec_tm)
        self.assertEqual(prec, true_prec)

class UnitTestEnumerate(unittest.TestCase):
    def test_enumerator(self):
        els = ['Cd','O']
        entries = GibbsSet(els,solids_data).entries
        reactions = EnumerateRxns(els,solids_data).rxns
        # number of reactions expected
        self.assertEqual(len(list(reactions)), 5)
        # all precursors/products are in allowed entries
        self.assertTrue(all(r in entries for r in reactions.entries))

class UnitTestTempEnv(unittest.TestCase):

    def test_mu(self):
        env_entries = [{'environment':'air', 'open':True},
                       {'environment':'inert', 'open':True}]           
        temperature = 300
        for entry in env_entries:
            if entry['environment'] not in ['air','inert']:
                self.assertEqual(TempEnvCorrections(entry['environment'], temperature).mu, None)
            elif entry['environment'] == 'air':
                self.assertAlmostEqual(TempEnvCorrections(temperature = temperature, environment = 'air').mu, 0.-0.040434717, places = 3)
            else:
                self.assertAlmostEqual(TempEnvCorrections(temperature = temperature, environment='inert').mu, -0.35749236, places=3)

    
    def test_entries_at_new_temp(self):
        reaction_set= EnumerateRxns(['Cd','O'],solids_data, temperature=300).rxns
        rxns_new = RxnsAtNewTempEnv(reaction_set, els = ['Cd','O'], solids_data = solids_data, new_temperature=1000, environment='air').corrected_reactions_at_temp()
        entries_new = rxns_new.entries
        new_temp_rxns = EnumerateRxns(['Cd','O'],solids_data, temperature=1000).rxns
        new_temp_rxns_corr = TempEnvCorrections(temperature=1000, environment='air').rxns_with_temp_env_correction(new_temp_rxns)
        new_temp_entries = new_temp_rxns_corr.entries
        for entry in entries_new:
            formula = entry.reduced_formula
            new_temp_entry = [entry for entry in new_temp_entries if entry.reduced_formula == formula][0]
            self.assertAlmostEqual(entry.energy, new_temp_entry.energy, places=3)


    def test_corrected_reactions_at_temp(self):
        reaction_set= EnumerateRxns(['Hf','O'],solids_data, temperature=300).rxns
        reaction_set_corrected = RxnsAtNewTempEnv(reaction_set, els = ['Hf','O'], solids_data = solids_data, new_temperature=1000, environment='air').corrected_reactions_at_temp()
        te = TempEnvCorrections(temperature = 1000, environment = 'air')
        reaction_list = reaction_set.get_rxns()
        for rxn in reaction_list:
            r_new = [entry for entry in reaction_set_corrected.get_rxns() if str(entry) == str(rxn)][0]
            corrections = []
            n_O = 0
            entries, coefficients = rxn.entries, rxn.coefficients
            for entry in entries:
                formula = entry.reduced_formula
                coefficient = coefficients[entries.index(entry)]
                if coefficient > 0:
                    if 'O' in formula:
                        n_O += CompTools(formula).amts['O']*coefficient        
                if te.environment_correction(formula) is not None:
                    corrections.append(te.environment_correction(formula).value * coefficient)
            reaction_energy += sum(corrections)
            if te.mu is not None:
                reaction_energy -= te.mu*n_O
            self.assertAlmostEqual(r_new.energy, reaction_energy, places=3)
                    
class UnitTestAnalyzeReactionSet(unittest.TestCase):
    def test_target_rxns(self):
        target = 'Hf1O2'
        environments = ['air','inert']
        for env in environments:
            rxns = EnumerateRxns(['Hf','O'],solids_data, temperature=300).rxns
            target_rxns = AnalyzeReactionSet(rxns, target = target, temperature = 300, environment= env).target_rxns()
            if env not in GAS_PARTIAL_PRESSURES:
                allowed_products = [target]
            else:
                allowed_products = [target] + [CompTools(p).clean for p in list(GAS_PARTIAL_PRESSURES[env].keys())]
            for rxn in target_rxns.get_rxns():
                self.assertTrue(CompTools(target).clean in rxn.products)
                self.assertTrue(all(CompTools(p).clean in allowed_products for p in rxn.products))
                if env not in GAS_PARTIAL_PRESSURES:
                    self.assertTrue(len(rxn.products)<=1)
                else:
                    self.assertTrue(len(rxn.products)<=2)

    def test_metrics_at_temp_env(self):
        target = 'Cd1S1'
        temperature = 1000
        rxns = EnumerateRxns(['Cd','S'],solids_data, temperature=temperature).rxns
        target_rxns = AnalyzeReactionSet(rxns, target = target, temperature = temperature, environment = 'air').target_rxns()
        ar = AnalyzeReactionSet(rxns, target = target, temperature = temperature, environment = 'air')
        metrics = ar.metrics_at_temp_env()
        for r in metrics:
            tr = [t for t in target_rxns if str(t) == r['rxn']][0]
            reactants, products = str.split(r['rxn'], ' -> ')
            products = [CompTools(p).clean for p in [products]]
            rxn_energy_hard = 0
            n_atoms = 0
            precursors = []
            for p in list(tr.entries):
                if CompTools(p.reduced_formula).clean in products:
                    rxn_energy_hard += p.energy
                    n_atoms += CompTools(p.reduced_formula).n_atoms
                else:
                    precursors.append(CompTools(p.reduced_formula).clean)
                    rxn_energy_hard -= p.energy
            rxn_energy_hard = rxn_energy_hard/n_atoms
            self.assertAlmostEqual(rxn_energy_hard, r['energy'], places= 3)
            prec_rxns = rxns.get_rxns_by_reactants(precursors)
            prec_rxns = [b for b in prec_rxns if str(b) != r['rxn']]
            competing_prec_rxn =  min(prec_rxns, key=lambda x: x.energy_per_atom)
            c1 = rxn_energy_hard - (competing_prec_rxn.energy_per_atom)
            self.assertAlmostEqual(c1, r['c1'], places = 3)
            precursors.append(target)
            secondary_rxns = list(rxns.get_rxns_by_reactants(precursors))
            # secondary reactions including target in precursors
            secondary_rxns = [r for r in secondary_rxns if CompTools(target).clean in [CompTools(p.reduced_formula).clean for p in r.reactants]]
            energies = []
            for sec in secondary_rxns:
                if sec.energy_per_atom < 0:
                    energies.append(sec.energy_per_atom)
            c2 = (sum(energies)/len(energies))*(-1)
            self.assertAlmostEqual(c2, r['c2'], places = 3)

class UnitTestRxnsAtNewTempEnv(unittest.TestCase):
    def test_corrected_reactions_at_temp(self):
        """
        Demonstrate the RxnsAtNewTempEnv to modify enumerated reactions from one temperature to a new temperatures 
        yields equivalent reactions to those generated using EnumerateRxns for a given temperature. 
        """
        rxns_at_300 = EnumerateRxns(['Cd','S'],solids_data, temperature=300).rxns
        rxns_changed_temp = RxnsAtNewTempEnv(rxns_at_300, ['Cd', 'S'], new_temperature= 1073, solids_data= solids_data).corrected_reactions_at_temp()
        rxns_at_1073 = EnumerateRxns(['Cd','S'],solids_data, temperature=1073).rxns
        rxns_at_1073 = TempEnvCorrections(temperature=1073).rxns_with_temp_env_correction(rxns_at_1073)
        self.assertEqual(len(list(rxns_changed_temp)), len(list(rxns_at_1073)))
        metrics_changed = AnalyzeReactionSet(reactions = rxns_changed_temp,
                                             target = 'Cd2S1', 
                                             temperature = 1073).metrics_at_temp_env()
        metrics_1073 = AnalyzeReactionSet(reactions = rxns_at_1073,
                                          target = 'Cd2S1',
                                          temperature= 1073).metrics_at_temp_env()
        for metric in metrics_changed:
            reaction_dict = get_reaction_dict_from_string(metric['rxn'])
            for entry in metrics_1073:
                print(entry)
                reaction_dict_1073 = get_reaction_dict_from_string(entry['rxn'])
                if set(reaction_dict_1073['reactants'])== set(reaction_dict['reactants']):
                    if set(reaction_dict_1073['products']) == set(reaction_dict['products']):
                        self.assertAlmostEqual(metric['energy'], entry['energy'], delta = 0.0001)
                        self.assertAlmostEqual(metric['c1'], entry['c1'], delta = 0.0001)
                        self.assertAlmostEqual(metric['c2'], entry['c2'], delta = 0.0001)
                        self.assertAlmostEqual(metric['gamma'], entry['gamma'], delta = 0.0001)


if __name__ == "__main__":
    unittest.main()