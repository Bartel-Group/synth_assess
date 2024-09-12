import os
from jobflow import run_locally
from rxn_network.flows.core import SynthesisPlanningFlowMaker
from rxn_network.jobs.core import GetEntrySetMaker
from rxn_network.jobs.core import CalculateCompetitionMaker
import numpy as np
from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json, write_json
from pymatgen.core.periodic_table import Element
from solidstatesynth.utils import *
# from utils import *



def get_chemical_potential(env, temp):
    std_pres = 1E5
    boltzmann_const = 8.617E-5
    # element = CompTools(element).clean
    # if element not in ['O', 'H', 'C1O2', 'H2O1','O2']:
    #     return 0
    # if element == 'O2':
    #     element = 'O'

    if env in ['nitrogen', 'argon', 'hydrogen', 'N2','Ar','H2', 'inert', 
                 'carbon monoxide', 'CO','carbon','platinum',]:
        pp_dict = {'O': 1E-6,
               'C1O2':1E-6,
               'H2O1':1E-6}
    elif env in ['oxygen','O2']:
        pp_dict = {'O': std_pres,
               'C1O2':1E-6,
               'H2O1':1E-6
               }
    else:
        pp_dict = {'O': 2E4,
               'C1O2':40,
               'H2O1':2.3E3
               }
    pres = pp_dict['O']
    mu = boltzmann_const * temp * np.log(pres / std_pres)
    return mu
    


def setup_flow(target, precursors=None, temp=300, env = 'air',entries=None, open = True):
    """
    This just defines the "FlowMaker" (what sequence of jobs will be run)

    In this case the flow is:
        1) collect and process entries
        2) enumerate reactions
        3) calculate competition
    """
    gesm = GetEntrySetMaker(temperature=temp, e_above_hull=0.050, filter_at_temperature=300)
    ccm = CalculateCompetitionMaker(calculate_chempot_distances=False)
    mu = get_chemical_potential(env,temp)
    # mu = 0
    if precursors == None or precursors == []:
        kwargs = {"exclusive_targets": True}
    else:
        if open:
            if 'O2' not in precursors:
                precursors.append('O2')
        kwargs = {"precursors": precursors, "exclusive_precursors": True, "exclusive_targets": True}

    extra_elements = set()
    if precursors != None and precursors != []:
        target_els = set(CompTools(target).els)
        prec_els = set()
        prec_els.update(*[CompTools(p).els for p in precursors])
        extra_elements = prec_els - target_els
    if open:
        flow_maker = SynthesisPlanningFlowMaker(get_entry_set_maker=gesm, basic_enumerator_kwargs=kwargs, minimize_enumerator_kwargs=kwargs,calculate_competition_maker=ccm, open_elem=Element('O'), chempots=[mu])
    else:
        flow_maker = SynthesisPlanningFlowMaker(get_entry_set_maker=gesm, basic_enumerator_kwargs=kwargs, minimize_enumerator_kwargs=kwargs,calculate_competition_maker=ccm)
    flow = flow_maker.make(target, entries=entries, added_elems=extra_elements)
    return flow


def parse_flow(flow):
    """
    This helps figure out what job ID is associated with each job in the flow
    """
    names_and_ids = {}
    for i in range(len(flow)):
        name = flow[i].name
        if "process" in name:
            name = "process"
        elif "enumerate" in name:
            name = "enumerate"
        elif "competition" in name:
            name = "competition"
        names_and_ids[name] = flow[i].uuid
    return names_and_ids


def parse_process_response(output):
    """
    in case you want to parse the first job (entry processing)
    """
    d = output.dict()
    entries = d["entries"].as_dict()
    d["entries"] = entries.copy()

    return d


def parse_enumerate_response(output):
    """
    in case you want to parse the second job (enumeration)
    """
    d = output.dict()
    rxns = d["rxns"]
    d["rxns"] = rxns.as_dict()
    d["enumerators"] = [enum.as_dict() for enum in d["enumerators"]]

    return d


def parse_competition_response(output):
    """
    the main parsing --> get the competition associated with each enumerated reaction

    Returns:
        {'rxns' : [{'rxn': str, 'c1': float, 'c2': float, 'entries' : [dict]}]}
            for every enumerated reaction
    """
    d = output.dict()
    rxns = d["rxns"]

    output = []
    rxns = rxns.get_rxns()
    for r in rxns:
        e = r.energy_per_atom
        r_str = str(r)
        r = r.as_dict()
        data = r["data"]
        # print(data)
        if "primary_competition" in data:
            c1 = data["primary_competition"]
        else:
            c1 = None
        if "secondary_competition" in data:
            c2 = data["secondary_competition"]
        else:
            c2 = None
        if c1!=None and c2!=None and e:
            gamma = 0.45 * c1 + 0.45 * c2 + 0.1*e
        else:
            gamma = None

        entries = r["entries"]
        output.append({"rxn": r_str, "c1": c1, "c2": c2, "entries": entries, "energy": e, "gamma": gamma})

    return {"rxns": output}


def run_flow(flow, fjson):
    """
    This is the computationally intensive bit so we'll write results to a json
    """
    # if os.path.exists(fjson) and not remake:
    #     return read_json(fjson)
    responses = run_locally(flow)

    output = {}

    names = parse_flow(flow)

    process_id = names["process"]
    process_response = responses[process_id][1].output
    output['process_entries'] = parse_process_response(process_response)
    # from pymatgen.core.composition import Composition
    # entries = output['process_entries']['entries']['entries']
    #print(entries)
    # formulas_in_my_entries = [Composition(e['composition']).reduced_formula for e in entries]
    # print(formulas_in_my_entries)
    # output['process_entries'] = parse_process_response(process_response)

    enumerate_id = names["enumerate"]
    enumerate_response = responses[enumerate_id][1].output
    # output['enumerators'] = parse_enumerate_response(enumerate_response)

    competition_id = names["competition"]
    # print('competition_id',competition_id)
    competitition_response = responses[competition_id][1].output
    output = parse_competition_response(competitition_response)

    return write_json(output, fjson)
    # return write_json(output, fjson)


def check_results(output):
    for r in output["rxns"]:
        print("\n")
        print(r["rxn"])
        print("c1 = ", r["c1"])
        print("c2 = ", r["c2"])

def get_competition_data(target,precursors,temperature, environment, open = True):
    fjson = "/Users/schle759/Mydrive/phd/research/flow_output.json"
    entries = None  # you can specify entries or let them be queried
    # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
    flow = setup_flow(target, precursors, temp = temperature, env = environment,entries=entries, open=open)
    names = parse_flow(flow)
    output = run_flow(flow, fjson=fjson)
    keys_to_get = ['rxn','c1','c2','energy','gamma']
    reactions = []
    for r in output["rxns"]:
        reactions.append({k: r[k] for k in keys_to_get})
        # print('reactions',reactions)
    # print(len(reactions))
    return reactions

def find_target_reaction(precursors, comp_data):
    if len(comp_data) == 0:
        # print('no target-reactions')
        return None
    for r in comp_data:
        rxn_dict = get_reaction_dict_from_string(r['rxn'])
        precursors = [CompTools(precursor).clean for precursor in precursors]
        r_prec = [CompTools(precursor).clean for precursor in rxn_dict['reactants']]
        # print('r_prec',r_prec)
        if set(r_prec) == set(precursors):
            return r
    # print('no prec-reactions')
    return None

def main():

    fjson = "/Users/schle759/Mydrive/phd/research/flow_output.json"

    entries = None  # you can specify entries or let them be queried
    # note: to query, you need to add your API_KEY to ~/.pmgrc.yaml (PMG_MAPI_KEY: < your API key >)
    # target = "La2Te1O6"
    precursors = ['O1', 'Te1', 'La2O3', 'O2Te1', 'La1','O3Te1','O1Te2']
    # , 'O3Te1', 'La2O5', 'O1Te2'
    flow = setup_flow(target="La2Te1O6", precursors = precursors, temp = 1300, entries=entries)
    names = parse_flow(flow)
    output = run_flow(flow, fjson=fjson)
    # competitions = get_competition_data(target,precursors,temperature = 1300)

    return output

if __name__ == "__main__":
    output = main()
