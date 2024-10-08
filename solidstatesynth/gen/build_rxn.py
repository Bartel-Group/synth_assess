import os
from jobflow import run_locally
# from rxn_network.flows.core import SynthesisPlanningFlowMaker
# from rxn_network.jobs.core import GetEntrySetMaker
# from rxn_network.jobs.core import CalculateCompetitionMaker
import numpy as np
from pydmclab.core.comp import CompTools
from solidstatesynth.extract.tm import get_updated_textmined_data
from pydmclab.utils.handy import read_json, write_json
from pymatgen.core.periodic_table import Element
from solidstatesynth.analyze.rxn import AnalyzeSynthesisRecipe
from solidstatesynth.analyze.compound import AnalyzeChemsys
from solidstatesynth.gen.metrics_calculation import EnumerateRxns, TargetRxns

# from solidstatesynth.analyze.compound import AnalyzeTarget
# from utils import *
DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"

class BuildRxn():
    def __init__(self, target, temperature=300, environment='air', with_theoretical=True):
        self.target = target
        self.temperature = temperature
        self.environment = environment
        self.mp_experimental = read_json(os.path.join(DATADIR, '241002_mp_experimental.json'))['data']
        self.stability_filter = 0.05
        self.with_theoretical = with_theoretical
        self.tm = get_updated_textmined_data()
        if with_theoretical:
            self.mp_data = read_json(os.path.join(DATADIR, '241002_mp_gd.json'))['data']
        else:
            self.mp_data = self.mp_experimental
        


    def get_precursors(self):
        """
        Returns a list of precursors appropriate for the desired target. First this function attempts to extract
        precursors from the text-mined dataset that are also in the materials project dataset. If the precursors
        listed cannot balance the reaction to make the target, the function will instead return all relevant 
        experimental precursors in the Materials Project dataset 

        """
        target_els = CompTools(self.target).els
        precursors = AnalyzeChemsys(target_els).possible_precursors(restrict_to_tm = True)
        balanceable = AnalyzeSynthesisRecipe(precursors=precursors,target=self.target, temperature=self.temperature, environment = self.environment).balanceable
        if not balanceable:
            precursors = AnalyzeChemsys(target_els).possible_precursors(restrict_to_tm = False)
        return precursors

    def build_target_rxns(self, found_precursors = None):
        """
        Args: option to specify precursors (as a list of strings)
        Returns a list of dictionaries with keys 'rxn', 'precursors', 'target', 'temperature', 'atmosphere', 'doi', 'mp'
        """
        # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
        if not found_precursors:
            found_precursors = self.get_precursors()
        print('precursors determined')
        target = self.target
        target_els = CompTools(target).els
        temperature = self.temperature
        environment = self.environment
        stability_filter = self.stability_filter
        with_theoretical = self.with_theoretical
        rxns = EnumerateRxns(els = target_els, temperature = temperature, with_theoretical = with_theoretical, stability_filter=stability_filter).rxns
        rxns_and_metrics =TargetRxns(precursors = found_precursors, target = target,reactions = rxns, temperature = temperature, environment = environment).metrics_at_temp_env()
        print('metrics calculated')
        return rxns_and_metrics


class AnalyzeRxnString():
    """
    Class to analyze a reaction string and determine if it is useful for a given target.
    This class is designed to function in conjunction with the reaction strings enumerated by
    the BuildRxn class.
    """
    def __init__(self, rxn_string):
        self.rxn_string = rxn_string

    @property
    def rxn_dict(self):
        """
        Args: reaction string
        Returns: dictionary with keys 'reactants' and 'products' where
        the values are lists of the reactants and products in the reaction
        *** IF CLEANABLE -- otherwise returns None ***
        Uses: this is the most usable reaction format from which to calculate
        dG_rxn using the pydmclab ReactionEnergy class-- get_dGrxn_at_T takes
        a reaction dictionary as an argument
        """
        rxn_string = self.rxn_string
        rxn_dict = {'precursors':[], 'products':[]}
        if '->' in rxn_string:
            sides = list(rxn_string.split(" -> "))
        elif '==' in rxn_string:
            sides = list(rxn_string.split(" == "))
            for i in range(2): 
                x = sides[i]
                entry_list = []   
                if "+" in x:
                    x = x.split(" + ")
                else:
                    x = [x]
                for entry in x:
                    if " " in entry:
                        if len(entry.split(" ")) != 1:
                            coefficient, entry = entry.split(" ")
                            if float(coefficient) > 0:
                                try:
                                    CompTools(entry).clean
                                    entry_list.append(entry)
                                except:
                                    return None
                    else:
                        entry_list.append(entry)
                if i == 0:
                    rxn_dict['precursors'] = entry_list
                else:
                    rxn_dict['products'] = entry_list
        return rxn_dict

    @property
    def products(self):
        return self.reaction_dict['products']

    def precursors(self):
        return self.reaction_dict['precursors']

    @property
    def has_gaseous_byproduct(self):
        """
        Returns:
            True if any of the products are gases
        """
        products = self.products
        for p in products:
            if AnalyzeCompound(p).is_gas:
                return True
        return False

    @property
    def has_gaseous_precursor(self):
        """
        Returns:
            True if any of the precursors are gases
        """
        precursors = self.precursors
        for p in precursors:
            if AnalyzeCompound(p).is_gas:
                return True
        return False


    def has_specified_precursors(self, desired_precursors):
        precursors = self.precursors
        desired_precursors = set([CompTools(p).clean for p in precursors])
        if all(p in desired_precursors for p in precursors):
            return True
        return False

    def is_useful_reaction(self, desired_target):
        """
        filters reactions based on requirement of no solid byproducts
        """
        products = self.products
        # use rxn string for specific precursors rather than a general list and to account for byproducts
        if desired_target in products:
            if len(products) <3:
                if len(products) < 2 or self.has_gaseous_byproduct:
                    return True  
        return False


class AnalyzeEnumeratedRxns():
    """
    Class to filter and optimize reactions for a given target generated by the BuildRxn class
    Takes in a list of reaction dictionaries with keys 'rxn', 'c1','c2','gamma', 'energy
    """
    def __init__(self, enumerated_rxns):
        self.enumerated_rxns = enumerated_rxns

    @property
    def rxn_list(self):
        return [r['rxn'] for r in self.enumerated_rxns]

    def filtered_rxns(self, desired_target):
        """
        Filters the reactions for those that are useful for the target as defined in the AnalyzeRxnString class
        """
        rxns = self.rxn_list
        filtered_rxns = []
        for r in rxns:
            if AnalyzeRxnString(r).is_useful_reaction(desired_target):
                filtered_rxns.append(r)
        return filtered_rxns
    
    
    def optimum_rxn(self, filter = 'gamma'):
        """
        Filter is defined as the key in the reaction dictionary that you want to optimize for
        identifying the best reaction. The default optimization filter is gamma.
        """
        rxns = self.enumerated_rxns
        if rxns:
            optimum_rxn = rxns[0]
            for r in rxns:
                if r[filter] < optimum_rxn[filter]:
                    optimum_rxn = r
        return optimum_rxn
        
        
                



def main():
    return 

if __name__ == "__main__":
    main()
