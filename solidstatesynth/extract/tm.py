""" 
The purpose of this module is to extract data from the text-mined database

Goal: map the text-mined dataset into simple dictionaries with easily accessible information

Don't worry about any thermo stuff here since this all gets taken care of in rxn_network (with or without our corrections)
"""

import os
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.query import MPRester

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"



def get_updated_textmined_data():
    """
    Returns: list of dictionaries of the form
    {"common":{"doi": ..., "temperature":..., "environment":...,
    "chemsys":..., "precursors":...}, "positive":{ "reaction": reaction,
    "target": target, "dG_rxn": dG_rxn}}
    with each item associated with an entry in "thermo" (a single reaction)
    """
    fjson = os.path.join(DATA_DIR, "240701_corrections.json")
    newdb = read_json(fjson)
    return newdb["data"]

def get_mp_cmpds():
    """
    Returns:
        list of compounds in MP
    """
    fjson = os.path.join(DATA_DIR, "241002_mp_ground_data.json")
    mp = read_json(fjson)['data']
    return list(mp.keys())


def get_tm_rxns(d, mp_cmpds, fjson=os.path.join("../data/tm_rxns.json"), remake=False):
    """
    Args:
        d (dict) : dictionary of text-mined data
            240626_no_corrections.json

        mp_cmpds (list) : list of compounds in MP

        fjson (str) : path to json file to save the output

        remake (bool) : if True, remake the json file

    Returns:
        [list of dictionaries where each dictionary is an observed reaction]

        for each reaction:
            {'rxn' : reaction string
            'precursors' : [list of precursor strings],
            'target': target string,
            'temperature' : int or None
            'atmosphere' : str or None,
            'doi' : str,
            'mp' : True if all precursors and target are in MP else False}
    """
    if not remake and os.path.exists(fjson):
        return read_json(fjson)["data"]

    out = {"data": []}
    for rxn in d:
        print(rxn)
        if ("positive" in rxn) and rxn["positive"]:
            precursors = rxn["common"]["precursors"]
            target = rxn["positive"]["target"]
            involved_cmpds = [p for p in precursors if precursors] + [target]
            involved_cmpds = list(set([c for c in involved_cmpds if c]))
            if all([c in mp_cmpds for c in involved_cmpds]):
                rxn["mp"] = True
            else:
                rxn["mp"] = False
            out["data"].append(
                {
                    "rxn": rxn["positive"]["reaction"],
                    "precursors": rxn["common"]["precursors"],
                    "target": rxn["positive"]["target"],
                    "temperature": rxn["common"]["temperature"],
                    "atmosphere": rxn["common"]["environment"],
                    "doi": rxn["common"]["doi"],
                    "mp": rxn["mp"],
                }
            )

    write_json(out, fjson)
    return out["data"]


def get_tm_precursors(
    tm_rxns,
    mp_cmpds,
    only_mp=True,
    fjson=os.path.join("../data/tm_precursors.json"),
    remake=False,
):
    """
    Args:
        tm_rxns (dict) : minimal dictionary of text-mined data
            generated with get_tm_rxns

        mp_cmpds (list) : list of compounds in MP

        only_mp (bool) : if True, only include precursors that are in MP

        fjson (str) : path to json file to save the output

        remake (bool) : if True, remake the json file


    Returns:
        list of precursors in the text-mined dataset that are also in MP (if only_mp=True)
    """
    if not remake and os.path.exists(fjson):
        return read_json(fjson)["data"]

    list_of_precursors = []
    for rxn in tm_rxns:
        if ("precursors" in rxn) and rxn["precursors"]:
            precursors = rxn["precursors"]
            list_of_precursors.extend(precursors)

    set_of_precursors = set(list_of_precursors)
    if only_mp:
        set_of_precursors = set_of_precursors.intersection(set(mp_cmpds))
    out = {"data": list(set_of_precursors)}
    write_json(out, fjson)
    return out["data"]


def get_tm_targets(
    tm_rxns,
    mp_cmpds,
    only_mp=True,
    fjson=os.path.join("../data/tm_targets.json"),
    remake=False,
):
    """
    Args:
        tm_rxns (dict) : minimal dictionary of text-mined data
            generated with get_tm_rxns

        mp_cmpds (list) : list of compounds in MP

        only_mp (bool) : if True, only include targets that are in MP

        fjson (str) : path to json file to save the output

        remake (bool) : if True, remake the json file

    Returns:
        list of targets in the text-mined dataset that are also in MP (if only_mp=True)
    """
    if not remake and os.path.exists(fjson):
        return read_json(fjson)["data"]

    list_of_targets = []
    for rxn in tm_rxns:
        if ("target" in rxn) and rxn["target"]:
            target = rxn["target"]
            list_of_targets.append(target)
    set_of_targets = set(list_of_targets)
    if only_mp:
        set_of_targets = set_of_targets.intersection(set(mp_cmpds))
    out = {"data": list(set_of_targets)}
    write_json(out, fjson)
    return out["data"]


def get_mp_computable_tm_rxns(
    tm_rxns,
):
    """
    Args:
        tm_rxns (dict) : minimal dictionary of text-mined data
            generated with get_tm_rxns

    Returns:
        list of reaction dictionaries in the text-mined dataset that are computable in MP
            this means all precursors and target are in MP
    """
    return [rxn for rxn in tm_rxns if rxn["mp"]]


def check():
    d = get_updated_textmined_data()
    print("%i entries in the text-mined dataset" % len(d))
    mp_cmpds = get_mp_cmpds()
    print("%i compounds in MP" % len(mp_cmpds))
    tm_rxns = get_tm_rxns(d, mp_cmpds)
    print("%i entries in the text-mined dataset after minimalizing" % len(tm_rxns))
    tm_rxns_computable_with_mp = get_mp_computable_tm_rxns(tm_rxns)
    print(
        "%i entries in the text-mined dataset computable with MP"
        % len(tm_rxns_computable_with_mp)
    )
    tm_precursors = get_tm_precursors(tm_rxns, mp_cmpds)
    print("%i precursors in the text-mined dataset are in MP" % len(tm_precursors))
    tm_targets = get_tm_targets(tm_rxns, mp_cmpds)
    print("%i targets in the text-mined dataset are in MP" % len(tm_targets))


def main():
    d = get_updated_textmined_data()
    mp_cmpds = get_mp_cmpds()
    tm_rxns = get_tm_rxns(d, mp_cmpds, remake=False)
    tm_precursors = get_tm_precursors(tm_rxns, mp_cmpds, remake=False)
    tm_targets = get_tm_targets(tm_rxns, mp_cmpds, remake=False)
    tm_rxns_in_mp = get_mp_computable_tm_rxns(tm_rxns)
    check()
    return d, tm_rxns, tm_precursors, tm_targets, tm_rxns_in_mp


if __name__ == "__main__":
    d, tm_rxns, tm_precursors, tm_targets, tm_rxns_in_mp = main()
