import os
from jobflow import run_locally
from rxn_network.flows.core import SynthesisPlanningFlowMaker
from rxn_network.jobs.core import GetEntrySetMaker
from rxn_network.jobs.core import CalculateCompetitionMaker
import numpy as np
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json, write_json
from pymatgen.core.periodic_table import Element
from solidstatesynth.analyze.rxn import AnalyzeRxn, AnalyzeRxnString
from solidstatesynth.analyze.compound import AnalyzeTarget
from solidstatesynth.gen.metrics_calculation import MetricsCalculator
from solidstatesynth.extract.mp import get_mp_experimental
# from solidstatesynth.analyze.compound import AnalyzeTarget
# from utils import *
DATADIR = "/Volumes/cems_bartel/projects/negative-examples/data"

class BuildRxn():
    def __init__(self, target, temp=300, env='air', with_theoretical=False):
        self.target = target
        self.temp = temp
        self.env = env
        self.with_theoretical = with_theoretical
        self.stability_filter = 0.1
        if with_theoretical:
            self.mp_data = read_json(os.path.join(DATADIR, '240925_mp_ground_data.json'))
        else:
            self.mp_data = read_json(os.path.join(DATADIR, '240926_mp_experimental.json'))


    def get_precursors(self):
        precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = True)
        balanceable = AnalyzeRxn(precursors=precursors,target=self.target, temperature=self.temp, atmosphere = self.env).balanceable
        if not balanceable:
            precursors = AnalyzeTarget(self.target).possible_precursors(restrict_to_tm = False, with_theoretical = self.with_theoretical)
        return precursors

    def build_target_rxns(self, precursors = None, open = True):
        # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
        if not precursors:
            precursors = self.get_precursors()
        print('precursors determined')
        target = self.target
        stability_filter = self.stability_filter
        with_theoretical = self.with_theoretical
        mc = MetricsCalculator(precursors=precursors,target = target, with_theoretical = with_theoretical, stability_filter=stability_filter)
        rxns_and_metrics =mc.metrics_at_temp_env(temp = self.temp, env = self.env)
        print('metrics calculated')
        return rxns_and_metrics
    

    def filtered_rxns(self, precursors = None, open = True):
        desired_target = self.target
        # print(target)
        reactions = self.build_target_rxns(precursors, open)
        print('reactions built')
        filtered_rxns = []
        for r in reactions:
            ars = AnalyzeRxnString(rxn=r['rxn'], atmosphere=self.env)
            if ars.is_useful_reaction(desired_target):
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
