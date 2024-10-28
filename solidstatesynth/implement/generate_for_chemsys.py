from solidstatesynth.gen.build_rxn import BuildRxn, AnalyzeRxnString, AnalyzeEnumeratedRxns
# from solidstatesynth.extract.tm import get_updated_textmined_data,get_mp_cmpds,get_tm_rxns
from pydmclab.core.comp import CompTools
from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeChemsys
from solidstatesynth.gen.metrics_calculation import EnumerateRxns, TargetRxns

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"

class ChemsysRxnGeneration():
    def __init__(self, chemsys,chemsys_targets):
        self.chemsys = chemsys
        self.chemsys_targets = targets
        self.chemsys_rxns = EnumerateRxns(chemsys).rxns


    def target_rxns_at_temp(self,target,temperature):
        precursors = BuildRxn(target = target).get_precursors()
        reactions = self.chemsys_rxns
        tr = TargetRxns(precursors = precursors, target = target, reactions = reactions)
        target_rxns = tr._find_target_rxns(temperature)
        new_temp_rxns = tr.get_rxns_at_new_temp_env(rxns=target_rxns,temp_new=temperature)
        return new_temp_rxns
    