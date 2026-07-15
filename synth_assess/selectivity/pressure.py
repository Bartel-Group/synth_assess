import os
from matcalc._eos import EOSCalc
from matcalc import load_fp
from pydmclab.core.query import MPQuery
from pydmclab.core.struc import StrucTools
from pydmclab.core.comp import CompTools
import numpy as np
import matplotlib.pyplot as plt

class StructureHofP:
    def __init__(
            self,
            struc,
            P_GPa,
            eos):
        
        self.struc = struc
        self.P_GPa = P_GPa
        self.eos = eos


    def eos_data(self):
        # EOS fitting with MatCalc (using MACE)
        eos = self.eos
        struc = self.struc
        try:
            data = eos.calc(struc)
        except:
            return [],[]
        V = np.array(data["eos"]["volumes"])
        E = np.array(data["eos"]["energies"])
        idx = np.argsort(V)
        V = V[idx]
        E = E[idx]
        if data['r2_score_bm']<0.95:
            print('low r2 score')
            return [], []
        return V,E
    
    def struc_V_H(self):

        GPA_TO_EVA3 = 1 / 160.21766208
        # Enthalpy at a given pressure determined given E(V)
        P = self.P_GPa * GPA_TO_EVA3  # convert once
        V,E = self.eos_data()
        if len(E) == 0:
            return None, None
        H_V = E + P*V
        # H_V = [E[i] + P * V[i] for i in range(len(E))]          # enthalpy at each sampled volume
        i_min = np.argmin(H_V)   # equilibrium volume
        return V[i_min], H_V[i_min]
    

class FormulaHofP:
    def __init__(
            self,
            formula,
            P_GPa,
            api = 'fDY65nUJRZpelqz70f36AHv3cxv2EpkY',
            eos = None,
            gs_struc_only = True):
        
        self.formula = formula
        self.P_GPa = P_GPa
        self.api = api
        self.gs_struc_only = gs_struc_only 
        if not eos:
            os.environ["HF_HUB_DISABLE_PROGRESS_BARS"] = "1"
            os.environ["TQDM_DISABLE"] = "1"
            model = load_fp("MACE")
            self.eos = EOSCalc(calculator=model, fmax=0.01)
        else:
            self.eos = eos

    def formula_strucs(self):
        mpq = MPQuery(self.api)
        if self.gs_struc_only:
            max_st = 1
        else:
            max_st = 10
        data = mpq.get_data(search_for = self.formula,
                        max_Ehull = 0.2,
                        max_polymorph_energy = 0.2,
                        max_strucs_per_cmpd = max_st,
                        max_sites_per_structure= 100,
                        additional_criteria = {"thermo_types":["GGA_GGA+U"], "theoretical":False})
        d_list = []
        for k in data:
            d_list.append({'structure':data[k]['structure'], 
                            'id':k,
                            'n_atoms': len(data[k]['structure']['sites'])})
        return d_list
    

    def formula_data(self):
        Hmin = 100
        Vmin = 0
        nsites = None
        struc_list = self.formula_strucs()
        if not struc_list:
            return None
        P_GPa = self.P_GPa
        eos = self.eos
        for entry in struc_list:
            st = StrucTools(entry['structure']).structure
            V,H = StructureHofP(st, P_GPa = P_GPa, eos = eos).struc_V_H()
            if not H:
                print('no st, continuing')
                continue
            if H < Hmin:
                Hmin = H
                Vmin = V
                nsites = len(st.sites)
        formula_entry = {'H_per_atom':Hmin/nsites, 'volume':Vmin, 'nsites':nsites}
        return formula_entry
        

