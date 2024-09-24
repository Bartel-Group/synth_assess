import os
from jobflow import run_locally
from rxn_network.flows.core import SynthesisPlanningFlowMaker
from rxn_network.jobs.core import GetEntrySetMaker
from rxn_network.jobs.core import CalculateCompetitionMaker
import numpy as np
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json, write_json
from pymatgen.core.periodic_table import Element
from solidstatesynth.gen.build_rxn_utils import setup_flow, parse_flow,run_flow
from solidstatesynth.analyze.rxn import AnalyzeRxn, AnalyzeRxnString
from solidstatesynth.analyze.compound import AnalyzeTarget
from solidstatesynth.gen.metrics_calculation import MetricsCalculator
# from solidstatesynth.analyze.compound import AnalyzeTarget
# from utils import *


class BuildRxn():
    def __init__(self, target, temp, env):
        self.target = target
        self.temp = temp
        self.env = env

    def get_precursors(self):
        precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = True)
        balanceable = AnalyzeRxn(precursors=precursors,target=self.target, temperature=self.temp, atmosphere = self.env).balanceable
        if not balanceable:
            precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = False)
        return precursors

    def build_target_rxns(self, precursors = None, open = True):
        fjson = "/Users/schle759/Mydrive/phd/research/flow_output.json"
        # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
        if not precursors:
            precursors = self.get_precursors()
        rxns_and_metrics = MetricsCalculator(precursors=precursors,targets = [self.target]).metrics_at_temp_env(temp = self.temp, env = self.env)

        # flow = setup_flow(self.target,temp = self.temp, env = self.env,open=open, precursors = precursors)
        # names = parse_flow(flow)
        # output = run_flow(flow, fjson=fjson)
        # keys_to_get = ['rxn','c1','c2','energy','gamma']
        # reactions = []
        # for r in output["rxns"]:
        #     reactions.append({k: r[k] for k in keys_to_get})
        return rxns_and_metrics
    

    def filtered_rxns(self, precursors = None, open = True):
        desired_target = self.target
        # print(target)
        reactions = self.build_target_rxns(precursors, open)
        filtered_rxns = []
        for r in reactions:
            ars = AnalyzeRxnString(rxn=r['rxn'], atmosphere=self.env)
            if ars.is_useful_reaction(desired_target):
                # if desired_target in ars.products:
                filtered_rxns.append(r)
        return filtered_rxns
    
    def optimum_rxn(self, filter = 'gamma',precursors = None, open = True):
        filtered_rxns = self.filtered_rxns(precursors, open)
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
