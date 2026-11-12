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
            P_GPa_list = [0,1,2,3,4,5,10],
            strain_list = [0.05,0.075,0.1,0.125]
            ):
        
        self.struc = struc
        self.P_GPa_list = P_GPa_list
        GPA_TO_EVA3 = 1 / 160.21766208
        self.P_list = [pressure*GPA_TO_EVA3 for pressure in P_GPa_list]
        self.strain_list = strain_list

    def eos_calc(self,max_abs_strain=0.05):
        os.environ["HF_HUB_DISABLE_PROGRESS_BARS"] = "1"
        os.environ["TQDM_DISABLE"] = "1"
        if self.model is None:
            self.model = load_fp("MACE")
        return EOSCalc(calculator=self.model, 
                       fmax=0.05, 
                       max_steps = 2000, 
                       max_abs_strain=max_abs_strain)



    def eos_data(self):
        # EOS fitting with MatCalc (using MACE)
        GPA_TO_EVA3 = 1 / 160.21766208
        struc = self.struc
        P = [pressure*GPA_TO_EVA3 for pressure in self.strain_list]
        for strain in [0.05,0.075,0.1,0.125]:
            eos = self.eos_calc(strain)
            try:
                data = eos.calc(struc)
            except:
                continue
            if data['r2_score_bm']<0.95:
                continue
            V = np.array(data["eos"]["volumes"])
            E = np.array(data["eos"]["energies"])
            idx = np.argsort(V)
            V, E = V[idx], E[idx]
            H = E[None, :] + P[:, None] * V[None, :]
            min_indices = np.argmin(H, axis=1)

            # Require minima to be away from boundaries
            if np.all(
                (min_indices > 1) &
                (min_indices < len(V) - 2)
            ):

                return V,E
        return [],[]
    
    def struc_V_H(self):

        # Enthalpy at a given pressure determined given E(V)
        P = np.array(self.P_list)
         # convert once

        V,E = self.eos_data()
        if len(E) == 0:
            return None,None
        H = E[None, :] + P[:, None] * V[None, :]

        min_indices = np.argmin(H, axis=1)

        V_min = V[min_indices]
        H_min = H[np.arange(len(P)), min_indices]

        if len(self.P_list) == 0:
            return V_min[0], H_min[0]

        return V_min, H_min

    

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
        elif self.P_GPa == 0:
            max_st = 1
        else:
            max_st = 10
        data = mpq.get_data(search_for = self.formula,
                        max_Ehull = 0.2,
                        max_polymorph_energy = 0.2,
                        max_strucs_per_cmpd = max_st,
                        max_sites_per_structure= 100,
                        properties= ['energy_per_atom'],
                        additional_criteria = {"thermo_types":["GGA_GGA+U"], "theoretical":False})
        d_list = []
        for k in data:
            d_list.append({'structure':data[k]['structure'], 
                            'id':k,
                            'n_atoms': len(data[k]['structure']['sites']),
                            'H_per_atom': data[k]['E_mp'],
                            'V': data[k]['volume']})
        return d_list
    

    def formula_data(self):
        Hmin = 100
        Vmin = 0
        nsites = None
        struc_list = self.formula_strucs()
        if not struc_list:
            return None
        P_GPa = self.P_GPa
        if P_GPa == 0:
            return {'H_per_atom':struc_list[0]['H_per_atom'], 
                    'volume':struc_list[0]['V'], 
                    'nsites':struc_list[0]['n_atoms']}
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
        

