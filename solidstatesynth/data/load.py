import os
from pydmclab.utils.handy import read_json


this_dir, this_filename = os.path.split(__file__)
DATA_PATH = os.path.join(this_dir, "data")

def mp_data():
    return read_json(os.path.join(DATA_PATH,'mp_solids_data.json'))

def tm_precursors():
    f = read_json(os.path.join(DATA_PATH,'tm_precursors.json'))
    return f['data']

def tm_rxns():
    f = read_json(os.path.join(DATA_PATH,'tm_rxns.json'))
    return f['data']

