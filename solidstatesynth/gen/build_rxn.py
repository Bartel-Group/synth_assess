import os
from jobflow import run_locally
from rxn_network.flows.core import SynthesisPlanningFlowMaker
from rxn_network.jobs.core import GetEntrySetMaker
from rxn_network.jobs.core import CalculateCompetitionMaker
import numpy as np
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json, write_json
from pymatgen.core.periodic_table import Element
from solidstatesynth.analyze.rxn import AnalyzeRxn, AnalyzeRxnDict
from solidstatesynth.analyze.compound import AnalyzeTarget
from solidstatesynth.gen.metrics_calculation import MetricsCalculator
from solidstatesynth.extract.mp import get_useful_mp_data
# from solidstatesynth.analyze.compound import AnalyzeTarget
# from utils import *
DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"

class BuildRxn():
    def __init__(self, target, temp=300, env='air', with_theoretical=True):
        self.target = target
        self.temp = temp
        self.env = env
        self.mp_experimental = read_json(os.path.join(DATADIR, '240926_mp_experimental.json'))['data']
        self.stability_filter = 0.05
        self.with_theoretical = with_theoretical
        if with_theoretical:
            self.mp_data = read_json(os.path.join(DATADIR, '240925_mp_ground_data.json'))['data']
        else:
            self.mp_data = self.mp_experimental
        self.useful_mp_data = get_useful_mp_data(self.mp_data)
        


    def get_precursors(self):
        """
        Returns a list of precursors appropriate for the desired target. First this function attempts to extract
        precursors from the text-mined dataset that are also in the materials project dataset. If the precursors
        listed cannot balance the reaction to make the target, the function will instead return all relevant 
        experimental precursors in the Materials Project dataset 

        """
        precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = True)
        balanceable = AnalyzeRxn(precursors=precursors,target=self.target, temperature=self.temp, atmosphere = self.env).balanceable
        if not balanceable:
            precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = False)
        return precursors

    def build_target_rxns(self, precursors = None):
        """
        Args: option to specify precursors (as a list of strings)
        Returns a list of dictionaries with keys 'rxn', 'precursors', 'target', 'temperature', 'atmosphere', 'doi', 'mp'
        """
        # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
        if not precursors:
            precursors = self.get_precursors()
        print('precursors determined')
        target = self.target
        temp = self.temp
        stability_filter = self.stability_filter
        with_theoretical = self.with_theoretical
        mc = MetricsCalculator(precursors=precursors,target = target, temperature = temp, with_theoretical = with_theoretical, stability_filter=stability_filter)
        rxns_and_metrics =mc.metrics_at_temp_env(env = self.env)
        print('metrics calculated')
        return rxns_and_metrics
    

    def get_reaction_dict_from_string(self, reaction_string):
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


    # Where to include get_reaction_dict_from_string

    def filtered_rxns(self, precursors = None):
        """
        Filters the reactions for those that are useful for the target as defined in the AnalyzeRxnString class
        """
        desired_target = self.target
        # print(target)
        reactions = self.build_target_rxns(precursors)
        print('reactions built')
        filtered_rxns = []
        for r in reactions:
            rxn_dict = self.get_reaction_dict_from_string(r['rxn'])
            ars = AnalyzeRxnDict(rxn_dict = rxn_dict)
            if ars.is_useful_reaction(desired_target):
                filtered_rxns.append(r)
        return filtered_rxns
    
    def optimum_rxn(self, filter = 'gamma',precursors = None):
        """
        Filter is defined as the key in the reaction dictionary that you want to optimize for
        identifying the best reaction. The default optimization filter is gamma.
        """
        filtered_rxns = self.filtered_rxns(precursors)
        optimized_rxn = None
        if filtered_rxns:
            optimized_rxn = filtered_rxns[0]
            for r in filtered_rxns:
                if r[filter] < optimized_rxn[filter]:
                    optimized_rxn = r
        return optimized_rxn
        
                



def main():
    return 

if __name__ == "__main__":
    main()
