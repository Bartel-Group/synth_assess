import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import CSS4_COLORS
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
from scipy.stats import iqr, wasserstein_distance, entropy
from scipy.optimize import minimize
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from pydmclab.plotting.utils import set_rc_params
from scipy.stats import binned_statistic
from mpl_toolkits.axes_grid1 import make_axes_locatable


set_rc_params()
COLORS = CSS4_COLORS # or some other palette
DATA_DIR = "/user/specifies/data/path"

class GammaDataProcessing:
    def __init__(
        self,
        tm_gamma_data,
        ):

        """
        class to determine optimum Gamma weights and return data with newly-computed/reweighted gamma.
        """
        self.gamma_data = tm_gamma_data
        self.tm_precursors = read_json(os.path.join(DATA_DIR,'tm_precursors.json'))
        ###### FIX THIS (STABILITY)
        self.stability = read_json(os.path.join(DATA_DIR, 'stability.json'))

    def data_reformatted(self):
        """
        helper function to format data for analysis. Redundant true reactions are aggregated in this step
        """
        data = self.gamma_data
        d_new_all = {key:{} for key in data}
        for key in data:
            for temperature in data[key]:
                if data[key][temperature]:
                    d_new = data[key][temperature][0]
                    if d_new:
                        d_new_all[key][temperature] = []
                        new_true = [e for e in d_new if e['true_rxn']]
                        if len(new_true) > 0:
                            for i in data[key][temperature]:
                                for j in range(len(i)):
                                    if i[j]['true_rxn'] == True:
                                        d_new[j]['true_rxn'] = True
                            for x in d_new:
                                x['temperature'] = int(float(temperature))
                        d_new_all[key][str(int(float(temperature)))] = d_new
        d_new_all= {key:d_new_all[key] for key in d_new_all if CompTools(key).n_els==3}
        return d_new_all

    def data_filtered(self, 
                      ternary = True):
        """
        helper function to filter data for analysis (removing non ternary oxides, empty entries)
        after data has been reformatted (refer to previous function). 
        """
        tm_precursors = self.tm_precursors
        data = self.data_reformatted()
        data_new = {key:{str(int(float(t))):[] for t in list(data[key].keys())} for key in data}
        for key in data:
            if ternary:
                if CompTools(key).n_els != 3:
                    continue
            if all([a not in CompTools(key).els for a in ['C', 'N', 'H']]):
                for temperature in data[key]:
                    for entry in data[key][temperature]:
                        rxn_dict = get_reaction_dict_from_string(entry['rxn'])
                        species = rxn_dict['reactants'] + rxn_dict['products']
                        species = [s for s in species if CompTools(s).clean != 'O2']
                        if all([7>CompTools(s).n_els for s in species]):
                            if len(rxn_dict['products'])>1:
                                entry['byproducts'] = True
                            else:
                                entry['byproducts'] = False
                            if tm_precursors:
                                if all([CompTools(i).clean in tm_precursors for i in rxn_dict['reactants']]):
                                    data_new[key][str(int(float(temperature)))].append(entry)

                            else:
                                data_new[key][str(int(float(temperature)))].append(entry)

            data_new[key] = {temperature: data_new[key][temperature] 
                                for temperature in data_new[key] if data_new[key][temperature]}

        data_new = {key: data_new[key] for key in data_new if data_new[key]}
        return data_new
    
    def data_list(self):
        """
        if minima is true, minimi
        """
        data = self.data_filtered()
        stability = self.stability

        data_entry_list = []
        for key in data:
            for temperature in data[key]:
                min_gamma = min([entry['gamma'] for entry in data[key][temperature]])
                min_entry = [entry for entry in data[key][temperature] if entry['gamma'] == min_gamma][0]
                min_entry['opt_rxn'] = True
                for entry in data[key][temperature]:
                    if temperature == 'null':
                        temperature = 1073.15
                    entry['temperature'] = int(float(temperature))
                    entry['target'] = CompTools(key).clean
                    entry['E_d'] = stability[CompTools(key).clean]
                    if 'opt_rxn' not in entry:
                        entry['opt_rxn'] = False
                    data_entry_list.append(entry)

        return data_entry_list

    
    
    def metric_data(self):
        data_list = self.data_list()
        metric_optima = []
        metric_true = []
        metric_all = []
        for entry in data_list:
            metric_all.append(entry['gamma'])
            if entry['true_rxn']:
                metric_true.append(entry['gamma'])
            if entry['opt_rxn']:
                metric_optima.append(entry['gamma'])
        return metric_optima, metric_true, metric_all

    def same_size_data(self):
        """
        Needed to compare true and optimum gamma when there are multiple true reactions associated
        with a given target and temperature (and optimum) 
        """
        data = self.data_filtered()
        gamma_optima = []
        gamma_true = []
        gamma_all = []
        for key in data:
            for temperature in data[key]:
                temp_data = data[key][temperature]
                if temp_data:
                    all_temp = [entry['gamma'] for entry in temp_data]
                    true = [entry['gamma'] for entry in temp_data if entry['true_rxn']]
                    opt = min(all_temp)
                    gamma_optima.extend([opt]*len(true))
                    gamma_true.extend(true)
        return gamma_optima, gamma_true, gamma_all

    def gamma_comp(self, weights):

        """
        helper function to identify optimum gamma weights. weights are input as a list [w1,w2,w3]
        """
        data_list = self.data_list()
        w1, w2, w3 = weights
        all_diffs = []
        
        for entry in data_list:
            true_rxn = [r for r in entry if r['true_rxn']]
            all_true_gamma = [w1 * r['c1'] + w2 * r['c2'] + w3 * r['energy'] for r in true_rxn]

            if all_true_gamma:
                all_true_gamma = [min(all_true_gamma)]
            else:
                continue
            gamma_opt = min(r['c1'] * w1 + r['c2'] * w2 + r['energy'] * w3 for r in entry)
            diff_list = [(t - gamma_opt)**2 for t in all_true_gamma]
            all_diffs.extend(diff_list)
        return np.mean(all_diffs)
    
    def optimize_gamma(self):
        constraints = ({'type': 'eq', 'fun': lambda w: np.sum(w) - 1})  # Constraint: w1 + w2 + w3 = 1
        bounds = [(0.0, 0.99)] * 3  # Ensure w1, w2, w3 are not 0 or 1
        data_list = self.data_list()
        # Define the grid of initial values (sum must be 1)
        grid_values = np.linspace(0.01, 0.99, 10)  # Adjust grid density as needed
        initial_guesses = [p for p in product(grid_values, repeat=3) if np.isclose(sum(p), 1.0, atol=1e-2)]
        best_result = 1e6
        best_w = None

        for guess in initial_guesses:
            result = minimize(self.gamma_comp, guess, args=(data_list,), bounds=bounds, constraints=constraints)
            if result.fun < best_result:
                best_result = result.fun
                best_w = result.x
                # print(best_w, best_result)

        return best_w  # Optimized w1, w2, w3
    
    def reweighted_gamma_data(self):
        data = self.data_filtered()
        w1, w2, w3 = self.optimize_gamma()
        for key in data:
            for temp in data[key]:
                for entry in data[key][temp]:
                    entry['gamma_new'] = w1*entry['c1'] + w2*entry['c2'] + w3*entry['energy']
        return data



def get_reaction_dict_from_string(reaction_string):
    """
    Args: reaction string
    Returns: dictionary with keys 'reactants' and 'products' where
    the values are lists of the reactants and products in the reaction
    *** IF CLEANABLE -- otherwise returns None ***
    Uses: this is the most usable reaction format from which to calculate
    dG_rxn using the pydmclab ReactionEnergy class-- get_dGrxn_at_T takes
    a reaction dictionary as an argument
    """
    reactant_list = []
    reactant_coeffs = []
    product_list = []
    product_coeffs = []
    if '->' in reaction_string:
        reactants, products = reaction_string.split(" -> ")
        # print(reactants)
    elif '==' in reaction_string:
        reactants, products = reaction_string.split(" == ")
    if "+" in reactants:
        reactants = reactants.split(" + ")
    else:
        reactants = [reactants]
    # print('reactants', reactants)
    for reactant in reactants:
        if " " in reactant:
            if len(reactant.split(" ")) != 1:
                coefficient, reactant = reactant.split(" ")
                if float(coefficient) > 0:
                    try:
                        CompTools(reactant).clean
                        reactant_list.append(reactant)
                        reactant_coeffs.append(coefficient)
                    except:
                        return None
        else:
            reactant_list.append(reactant)
            reactant_coeffs.append(1)
    if "+" in products:
        products = products.split(" + ")
    else:
        products = [products]
    for product in products:
        if " " in product:
            if len(product.split(" ")) != 1:
                coefficient, product = product.split(" ")
                if float(coefficient) > 0:
                    try:
                        CompTools(product).clean
                        product_list.append(product)
                        product_coeffs.append(coefficient)
                    except:
                        return None
        else:
            product_list.append(product)
            product_coeffs.append(1)
    # print(reactant_list,product_list)
    return {"reactants": reactant_list, "products": product_list, 
            "reactant_coeffs": reactant_coeffs, "product_coeffs": product_coeffs}



            
def main():
    solids_data = read_json(os.path.join(DATA_DIR, 'mp_solids_data.json'))
    tm_rxns = read_json(os.path.join(DATA_DIR, 'tm_data.json'))
    stability = read_json(os.path.join(DATA_DIR, 'stability.json'))
    # prediction data
    tm_gamma_data = read_json("path/to/metrics/file")
    # user specified?
    data_tm = GammaDataProcessing(tm_gamma_data).data_filtered()
    tm_list = GammaDataProcessing(tm_gamma_data).data_list()
    tm_true = [entry for entry in tm_list if entry['true_rxn']]
    tm_opt = [entry for entry in tm_list if entry['opt_rxn']]
    # pred_dict = synth_scores(synth_pred,sd, data_tm)
    # # write_json(pred_dict, os.path.join(DATA_DIR,'tm_predictions.json'))
    return solids_data, tm_rxns, stability, tm_gamma_data, data_tm, tm_list, tm_true, tm_opt

if __name__ == "__main__":
   solids_data, tm_rxns, stability, tm_gamma_data, tm_list, tm_true, tm_opt = main() 
