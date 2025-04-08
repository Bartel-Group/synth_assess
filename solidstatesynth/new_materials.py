from pydmclab.utils.handy import read_json
import os
from pydmclab.core.comp import CompTools
from solidstatesynth.gen.rxn_metrics import EnumerateRxns, RxnsAtNewTempEnv, AnalyzeReactionSet
DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"


def get_metrics(target, temperature, solids_data, gen_data = None, is_gen = False):
    if is_gen:
        gen_formula = target
    else:
        gen_formula = None
    r = EnumerateRxns(els = CompTools(target).els, temperature = 300, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula).rxns
    r = RxnsAtNewTempEnv(reaction_set = r, els = CompTools(target).els, new_temperature = temperature, solids_data = solids_data, gen_data = gen_data, gen_formula = gen_formula).corrected_reactions_at_temp()
    print(r.entries)
    t = AnalyzeReactionSet(target = target, reactions = r, temperature = temperature)
    return t.metrics_at_temp_env()  


def get_new_target_dict(target):
    temp_list = [600,900,1200,1500,1800]
    t_dict = [{'target':target, 'temperature':t} for t in temp_list]
    return t_dict

def get_new_solids(new_mat_list, metric_chemsystems):
    in_tm = []
    in_tm_dict = {}
    not_tm = []
    not_tm_dict = {}
    for mat in new_mat_list:
        els = CompTools(mat).els
        chemsys = CompTools(mat).chemsys
        if 'O' not in els:
            continue
        elif len(els) >4:
            continue
        elif len(els)<3:
            continue
        elif any([i in els for i in ['C','H']]):
            continue
        # elif CompTools(mat).clean in metric_targets:
        #     continue
        elif chemsys in metric_chemsystems:
            in_tm.append(mat)
        elif len(els) == 3:
            set_cs = set(chemsys.split('-'))
            for c in metric_chemsystems:
                if set_cs.issubset(set(c.split('-'))):
                    in_tm.append(mat)
            if mat not in in_tm:
                not_tm.append(mat)
        else:
            not_tm.append(mat)
    for entry in in_tm:
        t_dict = get_new_target_dict(entry)
        chemsys = CompTools(entry).chemsys
        chemsys = chemsys.replace('-',',')
        if chemsys not in in_tm_dict:
            in_tm_dict[chemsys]=[]
        in_tm_dict[chemsys].extend(t_dict)
    for entry in not_tm:
        t_dict = get_new_target_dict(entry)
        chemsys = CompTools(entry).chemsys
        chemsys = chemsys.replace('-',',')
        if chemsys not in not_tm_dict:
            not_tm_dict[chemsys]=[]
        not_tm_dict[chemsys].extend(t_dict)
    return in_tm_dict, not_tm_dict
    # return in_tm, not_tm


def main():
    solids_data = read_json(os.path.join(DATA_DIR, 'solids_data.json'))
    metrics = read_json(os.path.join(DATA_DIR, 'march_all_metrics.json'))
    metric_chemsystems = []
    for entry in metrics:
        cs = CompTools(list(entry.keys())[0]).chemsys
        if cs not in metric_chemsystems:
            metric_chemsystems.append(cs)
    gen_data = read_json(os.path.join(DATA_DIR, '250324_new_materials.json'))
    gen_data = {key: {CompTools(k).clean: gen_data[key][k] for k in gen_data[key]} for key in gen_data}
    gen_mat_filtered = read_json(os.path.join(DATA_DIR, 'gen_mat_new.json'))
    new_mat_list = list(set([k for key in gen_data for k in gen_data[key].keys()]))
    return solids_data, metric_chemsystems, gen_data, gen_mat_filtered,new_mat_list

if __name__ == "__main__":
    solids_data, metric_chemsystems, gen_data, gen_mat_filtered, new_mat_list = main()