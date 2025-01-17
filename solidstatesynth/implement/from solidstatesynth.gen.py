from solidstatesynth.gen.build_rxn import BuildRxn, AnalyzeRxnString, AnalyzeEnumeratedRxns
from solidstatesynth.extract.tm import get_updated_textmined_data,get_mp_cmpds,get_tm_rxns
from pydmclab.core.comp import CompTools
from solidstatesynth.analyze.compound import AnalyzeCompound, AnalyzeChemsys
from solidstatesynth.gen.metrics_calculation import EnumerateRxns

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"
common_chemsystems = ['Bi-Fe-O', 'Li-O-Ti', 'Li-Mn-O', 'Co-Li-O', 'O-Sr-Ti', 'Ba-O-Ti', 'Al-O-Y']

def tm_rxns():
    """
    Get all thermomechanical reactions from the Materials Project.
    """
    tm = get_updated_textmined_data()
    cmpds = get_mp_cmpds()
    tm_rxns = get_tm_rxns(tm, cmpds)
    return tm_rxns

def relevant_tm_rxns():
    """
    Get all thermomechanical reactions from the Materials Project that are relevant to the negative example generation.
    """
    rxns = tm_rxns()
    relevant_tm_rxns = []
    for rxn in rxns:
        target = rxn['target']
        target_els = CompTools(target).els
        if len(target_els)==3:
            if 'O' in target_els:
                environment = rxn['atmosphere']
                if environment == 'air' or environment == None:
                    relevant_tm_rxns.append(rxn)
    return relevant_tm_rxns

def chemsys_reactions(els):
    reactions = EnumerateRxns(els).rxns
    return reactions


def competitions_for_tm_entry(tm_entry):
    tm_target = tm_entry['target']
    tm_precursors = AnalyzeRxnString(tm_entry['rxn']).precursors
    tm_environment = tm_entry['atmosphere']
    if not tm_environment:
        tm_environment = 'air'
    tm_temperature = tm_entry['temperature']
    if not tm_temperature:
        tm_temperature = 1073
    tm_rxns = BuildRxn(target = tm_target, environment= tm_environment, temperature=tm_temperature).build_target_rxns()
    actual_rxn = None
    optimum = AnalyzeEnumeratedRxns(tm_rxns).optimum_rxn(tm_target)
    for rxn in tm_rxns:
        ar = AnalyzeRxnString(rxn['rxn'])
        if ar.has_specified_precursors(tm_precursors):
            actual_rxn = rxn
            return optimum, actual_rxn
    return optimum, actual_rxn

def get_competitions(rxns):
    competitions = []
    for tm_rxn in rxns:
        optimum, actual_rxn = competitions_for_tm_entry(tm_rxn)
        if actual_rxn:
            competitions.append({'optimum':optimum, 'actual':actual_rxn})
    return competitions

def main():
    tm = relevant_tm_rxns()
    return tm

if __name__ == "__main__":
    tm = main()