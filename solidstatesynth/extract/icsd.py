from pydmclab.core.comp import CompTools
from pydmclab.utils.handy import read_json, write_json

DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"



def get_new_target_dict(target):
    temp_list = [600,900,1200,1500,1800]
    t_dict = [{'target':target, 'temperature':t} for t in temp_list]
    return t_dict

def get_icsd_solids(solids, metrics):
    metric_targets = list(metrics.keys())
    metric_targets = [CompTools(t).clean for t in metric_targets]
    in_tm = []
    in_tm_dict = {}
    not_tm = []
    not_tm_dict = {}
    metric_chemsystems = list(set([CompTools(t).chemsys for t in metric_targets]))
    metric_chemsystems.extend(['Li-Mn-Ni-O','Co-Li-Mn-O','Li-O-P-V', 'Co-Li-O-P', 'Fe-Li-O-P',
                               'Li-Mn-O-P','Li-Mn-O-Ti','Fe-Li-O-Ti','Co-Li-O-Ti','Bi-Li-O-Ti',
                               'La-Li-Mn-O','Ca-Fe-O-P','Fe-Li-Mn-O','Li-Ni-O-V','Fe-Li-O-Si',
                               'Li-Ni-O-P','Fe-Na-O-P','B-Li-Mn-O','Cu-Li-O-V','Mn-Na-O-V',
                               'Li-Mn-O-Si','Cr-Li-O-Si','Cr-Li-O-Ti','Fe-O-P-Sr','Li-Mn-O-V',
                               'Fe-Li-O-V'])
    for key in solids:
        els = CompTools(key).els
        chemsys = CompTools(key).chemsys
        if 'O' not in els:
            continue
        elif len(els) >4:
            continue
        elif len(els)<3:
            continue
        elif any([i in els for i in ['C','H']]):
            continue
        elif CompTools(key).clean in metric_targets:
            continue
        elif chemsys in metric_chemsystems:
            in_tm.append(key)
        elif len(els) == 3:
            set_cs = set(chemsys.split('-'))
            for c in metric_chemsystems:
                if set_cs.issubset(set(c.split('-'))):
                    in_tm.append(key)
            if key not in in_tm:
                not_tm.append(key)
        else:
            not_tm.append(key)
    for entry in in_tm:
        t_dict = get_new_target_dict(entry)
        chemsys = CompTools(entry).chemsys
        chemsys = chemsys.replace('-',',')
        if chemsys not in in_tm_dict:
            in_tm_dict[chemsys]=[]
        in_tm_dict[chemsys].append(t_dict)
    for entry in not_tm:
        t_dict = get_new_target_dict(entry)
        chemsys = CompTools(entry).chemsys
        chemsys = chemsys.replace('-',',')
        if chemsys not in not_tm_dict:
            not_tm_dict[chemsys]=[]
        not_tm_dict[chemsys].append(t_dict)
    return in_tm_dict, not_tm_dict