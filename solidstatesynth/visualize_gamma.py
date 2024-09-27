from rxn_network.reactions.basic import BasicReaction
from pydmclab.utils.handy import read_json, write_json
import os
from pymatgen.core.composition import Composition
from solidstatesynth.core.utils import get_reaction_dict_from_string

gammas_min = {}

PATH_air = "/Volumes/cems_bartel/projects/negative-examples/data/simon/240922_calculations_air"
PATH_inert = "/Volumes/cems_bartel/projects/negative-examples/data/simon/240922_calculations_inert"
PATH_competition = "/Volumes/cems_bartel/projects/negative-examples/data"
#"Volumes\\cemsbartel\\\\DMC Lab\\ICSD\\240922_calculations_air"


def get_min_gamma_data(PATH):
    target_min_gamma = {}
    for filename in os.listdir(PATH):
        metrics_data_air = read_json(os.path.join(PATH, filename))
        # metrics_data_inert = read_json(os.path.join(PATH_inert, filename))
        for j in metrics_data_air:
            for i in metrics_data_air[j]:
                rxn = BasicReaction.from_string(i["rxn"])
                target = set(rxn.products) - {Composition("O2"), Composition("CO2")}
                target = target.pop()
                if target not in target_min_gamma:
                    target_min_gamma[target] = i
                    target_min_gamma[target]["temp"] = j
                elif target_min_gamma[target]["gamma"] >= i["gamma"]:
                    target_min_gamma[target] = i
                    target_min_gamma[target]["temp"] = j    
                gammas_min.update(target_min_gamma)
    return target_min_gamma

def get_gamma_data(PATH):
    target_gamma = {}
    for filename in os.listdir(PATH):
        metrics_data_air = read_json(os.path.join(PATH, filename))
        # metrics_data_inert = read_json(os.path.join(PATH_inert, filename))
        for j in metrics_data_air:
            for i in metrics_data_air[j]:
                rxn = BasicReaction.from_string(i["rxn"])
                target = set(rxn.products) - {Composition("O2"), Composition("CO2")}
                target = target.pop()
                if target not in target_gamma:
                    target_gamma[target] = i
                    target_gamma[target]["temp"] = [j]
                elif j not in target_gamma[target]["temp"]:
                    target_gamma[target]["temp"].append(j)
                gammas_min.update(target_gamma)
    return target_gamma


def get_competition_data():
    fjson = os.path.join(PATH_competition, "240705_mcdermott_data.json")
    competition = read_json(fjson)
    for entry in competition['Sheet1']:
        rxn_dict = get_reaction_dict_from_string(entry["reaction"])
        precursors = rxn_dict["reactants"]
        target = [entry for entry in rxn_dict["products"] if entry not in ["C1O2","CO2","O2","O1","H2O1","H2O"]]
        entry["precursors"] = precursors
        entry["target"] = target
    # competition = [entry for entry in competition['Sheet1'] if 'open' not in entry]
    return competition


def main():
    gammas_air = get_min_gamma_data(PATH_air)
    gammas_air_all = get_gamma_data(PATH_air)
    competition_data = get_competition_data()
    competition_targets = list(set([Composition(entry["target"][0]) for entry in competition_data['Sheet1']]))
    return gammas_air, gammas_air_all, competition_data, competition_targets


if __name__ == "__main__":
    (gammas_air,
     gammas_air_all, 
     competition_data,
     competition_targets,
     ) = main()