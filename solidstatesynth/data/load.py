import os
from pydmclab.utils.handy import read_json


this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "data")
results_DATA_PATH = os.path.join(this_dir, "results_data")

def mp_data():
    return read_json(os.path.join(DATA_PATH,'mp_solids_data.json'))

def tm_precursors():
    f = read_json(os.path.join(DATA_PATH,'tm_precursors.json'))
    return f['data']

def tm_rxns():
    f = read_json(os.path.join(DATA_PATH,'tm_rxns.json'))
    return f['data']

### results (use for plotting)
def tm_entries():
    f = read_json(os.path.join(results_DATA_PATH,'tm_entries.json'))
    return f['data']

def gen_data():
    f = read_json(os.path.join(results_DATA_PATH,'gen_mat_pred_gamma.json'))
    return f

def tm_rxns_with_gamma():
    f = read_json(os.path.join(results_DATA_PATH,'tm_rxns_with_gamma.json'))
    return f
