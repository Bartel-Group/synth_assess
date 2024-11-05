""" 
The purpose of this module is to retrieve MP data
"""

from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.query import MPQuery

import os


def get_mp_icsd_cmpds():
    return read_json("../data/mp_icsd_cmpds.json")["formulas"]


def get_all_mp_data(fjson="../data/240924_mp_icsd_data.json", remake=False):
    if os.path.exists(fjson) and not remake:
        return read_json(fjson)
    mpq = MPQuery()
    mp_icsd_formulas = get_mp_icsd_cmpds()
    data = mpq.get_data(
        search_for=mp_icsd_formulas,
        max_strucs_per_cmpd=1,
        max_sites_per_structure=None,
        max_Ehull=0.2,
        include_structure=False,
    )
    out = write_json(data, fjson)
    return out


def main():
    mp = get_all_mp_data()
    return mp


if __name__ == "__main__":
    main()
