import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import CSS4_COLORS
from pydmclab.utils.handy import read_json, write_json
from pydmclab.core.comp import CompTools
from scipy.stats import iqr, wasserstein_distance, entropy
from scipy.optimize import minimize, differential_evolution
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import product
from pydmclab.plotting.utils import set_rc_params
from pydmclab.core.energies import ChemPots, FormationEnergy
from pydmclab.core.hulls import GetHullInputData, AnalyzeHull


set_rc_params()
COLORS = CSS4_COLORS # or some other palette
DATA_DIR = "/Volumes/cems_bartel/projects/negative-examples/data"
DATA_DIR_metrics = "/Volumes/cems_bartel/projects/negative-examples/data/metrics"

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

def data_reformatted(data, is_icsd = False):
    d_new_all = {key:{} for key in data}
    for key in data:
        for temperature in data[key]:
            if data[key][temperature]:
                d_new = data[key][temperature][0]
                if d_new:
                    if not is_icsd:
                        d_new_all[key][temperature] = []
                        new_true = [e for e in d_new if e['true_rxn']]
                        if len(new_true) > 0:
                            for i in data[key][temperature]:
                                for j in range(len(i)):
                                    if i[j]['true_rxn'] == True:
                                        d_new[j]['true_rxn'] = True
                    d_new_all[key][temperature] = d_new
    return d_new_all


def data_filtered(data_r):
    data_new = {key:{t:[] for t in list(data_r[key].keys())} for key in data_r}
    for key in data_r:
        if all([a not in CompTools(key).els for a in ['C', 'N', 'H']]):
            for temperature in data_r[key]:
                for entry in data_r[key][temperature]:
                    rxn_dict = get_reaction_dict_from_string(entry['rxn'])
                    species = rxn_dict['reactants'] + rxn_dict['products']
                    species = [s for s in species if CompTools(s).clean != 'O1']
                    if all([5>CompTools(s).n_els > 1 for s in species]):
                        if len(rxn_dict['products'])>1:
                            entry['byproducts'] = True
                        else:
                            entry['byproducts'] = False
                        data_new[key][temperature].append(entry)
        data_new[key] = {temperature: data_new[key][temperature] 
                            for temperature in data_new[key] if data_new[key][temperature]}
    data_new = {key: data_new[key] for key in data_new if data_new[key]}
    return data_new

def data_gamma_new(data_rf, w1, w2, w3):
    for key in data_rf:
        for temp in data_rf[key]:
            for entry in data_rf[key][temp]:
                entry['gamma_new'] = w1*entry['c1'] + w2*entry['c2'] + w3*entry['energy']
    return data_rf

def get_ed_at_temp(target, temperature, solids_data):
    target = CompTools(target).clean
    chemsys = CompTools(target).chemsys
    t_els = chemsys.split('-')
    t_entries = list(set([k for k in solids_data if all([e in t_els for e in CompTools(k).els])]))
    energy_dict = {t: {'Ef':None} for t in t_entries}
    for t in t_entries:
        Ef = solids_data[t]['formation_energy_per_atom']
        vol = solids_data[t]['volume']/solids_data[t]['nsites']
        t_n = int(round(temperature, -2))
        if t_n >2000:
            t_n = 2000
        if t_n < 300:
            t_n = 300
        mus = ChemPots(temperature=t_n).chempots
        Ef_at_T = FormationEnergy(formula=t, Ef=Ef, chempots=mus, atomic_volume=vol).dGf(temperature=int(temperature))
        energy_dict[t]['Ef'] = Ef_at_T
    h = GetHullInputData(energy_dict, formation_energy_key='Ef').hullin_data(fjson = 'hull_t.json', remake=True)
    a = AnalyzeHull(hullin_data=h, chemical_space= chemsys.replace('-','_'))
    return a.decomp_energy(target)





def get_metric_data(metric, data_list, is_icsd = False, byproducts = None):
    metric_optima = []
    metric_true = []
    metric_all = []
    for entry in data_list:
        if not is_icsd:
            metric_all.append(entry[metric])
            if entry['true_rxn']:
                metric_true.append(entry[metric])
        if entry['opt_rxn']:
            metric_optima.append(entry[metric])
    return metric_optima, metric_true, metric_all

def gamma_comp(weights, entries):
    w1, w2, w3 = weights
    all_diffs = []
    
    for entry in entries:
        true_rxn = [r for r in entry if r['true_rxn']]
        all_true_gamma = [w1 * r['c1'] + w2 * r['c2'] + w3 * r['energy'] for r in true_rxn]
        
        gamma_opt = min(r['c1'] * w1 + r['c2'] * w2 + r['energy'] * w3 for r in entry)
        diff_list = [t - gamma_opt for t in all_true_gamma]
        all_diffs.extend(diff_list)
    
    return np.mean(all_diffs)

def optimize_gamma(entries):
    constraints = ({'type': 'eq', 'fun': lambda w: np.sum(w) - 1})  # Constraint: w1 + w2 + w3 = 1
    bounds = [(0.0, 0.99)] * 3  # Ensure w1, w2, w3 are not 0 or 1
    
    # Define the grid of initial values (sum must be 1)
    grid_values = np.linspace(0.1, 0.8, 4)  # Adjust grid density as needed
    initial_guesses = [p for p in product(grid_values, repeat=3) if np.isclose(sum(p), 1.0)]
    
    best_result = None
    best_w = None

    for guess in initial_guesses:
        result = minimize(gamma_comp, guess, args=(entries,), bounds=bounds, constraints=constraints)
        if best_result is None or result.fun < best_result.fun:
            best_result = result
            best_w = result.x

    return best_w  # Optimized w1, w2, w3

def get_formatted_data(file_name, optimize = True, is_icsd = False):
    data = read_json(os.path.join(DATA_DIR, file_name))
    data = {key: entry[key] for entry in data for key in entry}
    data = data_reformatted(data, is_icsd)
    data = data_filtered(data)
    if optimize:
        data_list = []
        for key in data:
            for temp in data[key]:
                data_list.append(data[key][temp])
        print(data_list[0])
        w1, w2, w3 = optimize_gamma(data_list)
        print('opt', w1,w2,w3)
    else:
        w1, w2, w3 = 0.093445, 0.63379, 0.272765
    print(w1,w2,w3)
    data_n = data_gamma_new(data, w1,w2,w3)
    return data_n

def data_entries(data, solids_data, minima = False, gamma = 'gamma_new', filename = 'data_entries.json', Ed = False):
    p = os.path.join(DATA_DIR,filename)
    if os.path.exists(p):
        return read_json(p)['data']
    data_entry_list = []
    for key in data:
        if minima == False:
            for temperature in data[key]:
                min_gamma = min([entry[gamma] for entry in data[key][temperature]])
                min_entry = [entry for entry in data[key][temperature] if entry[gamma] == min_gamma][0]
                min_entry['opt_rxn'] = True
                for entry in data[key][temperature]:
                    entry['temperature'] = float(temperature)
                    entry['target'] = key
                    if 'opt_rxn' not in entry:
                        entry['opt_rxn'] = False
                    if Ed:
                        try:
                            entry['E_d'] = get_ed_at_temp(key, float(temperature), solids_data=solids_data)
                        except:
                            print(key)
                            entry['E_d'] = None
                    data_entry_list.append(entry)

        else:
            mins = []
            for temperature in data[key]:
                min_gamma = min([entry[gamma] for entry in data[key][temperature]])
                min_entry = [entry for entry in data[key][temperature] if entry[gamma] == min_gamma][0]
                min_entry['temperature'] = float(temperature)
                min_entry['target'] = key
                mins.append(min_entry)
            min_min = min([entry[gamma] for entry in mins])
            min_entry = [entry for entry in mins if entry[gamma] == min_min]
            data_entry_list.append(min_entry)
    f = write_json({'data':data_entry_list},p)
    return read_json(p)['data']

def data_by_byproduct(data_entry_list, byproducts):
    d_list = []
    for entry in data_entry_list:
        rxn = entry['rxn']
        rxn_dict = get_reaction_dict_from_string(rxn)
        if byproducts == 'O2':
            if 'O2' in rxn_dict['products'] + rxn_dict['reactants']:
                d_list.append(entry)
        elif byproducts == None:
            if len(rxn_dict['products']) == 1:
                d_list.append(entry)
        elif byproducts == 'CO2':
            if 'C1O2' in rxn_dict['products'] or 'CO2' in rxn_dict['products']:
                d_list.append(entry)
        elif byproducts == 'H2O':
            if 'H2O1' in rxn_dict['products'] or 'H2O' in rxn_dict['products']:
                d_list.append(entry)
        elif byproducts == 'any':
            if len(rxn_dict['products'])>1:
                d_list.append(entry)
    return d_list



def get_metric_data_same_size(metric, data):
    metric_optima = []
    metric_true = []
    metric_all = []
    # min_data = []
    for key in data:
        for temperature in data[key]:
            temp_data = data[key][temperature]
            if temp_data:
                all_temp = [entry[metric] for entry in temp_data]
                true = [entry[metric] for entry in temp_data if entry['true_rxn']]
                opt = min(all_temp)
                metric_optima.extend([opt]*len(true))
                # min_data_o = [entry for entry in temp_data if entry[metric] == min(metric_all)]
                metric_true.extend(true)
    return metric_optima, metric_true, metric_all



# Example histogram data (replace with actual histograms)

def metric_histogram(ax, metric, data, is_icsd = False):
    if not is_icsd:
        metric_optima, metric_true, metric_all = get_metric_data(metric, data)
    else:
        metric_optima = get_metric_data(metric, data, is_icsd=True)
    # fig, ax1 = plt.subplots()
    ax.set_xlabel(f"{metric}")
    # ax.set_xlabel("$\Gamma$")
    sdev_all = np.std(metric_all)
    if metric == 'c1':
        ax.set_xlim(-0.5,0.5)
        # ax.set_xlim(np.mean(metric_all) - 1.5*sdev_all, np.mean(metric_all) + 2*sdev_all)
        ax.set_ylim(0,10)
        ymax = 15
        b_opt = 60
        b_true = 80
        b_all = 170
        zorder_big = 2
        zorder_small = 1
    if metric == 'c2':
        ax.set_xlim(0,0.7)
        ax.set_ylim(0,12)
        ymax = 20
        b_opt = 180
        b_true = 350
        b_all = 800
        zorder_big = 2
        zorder_small = 1
    if metric == 'energy':
        ax.set_xlim(-0.5, 0.3)
        ax.set_ylim(0,9)
        ymax = 14
        b_opt = 200
        b_true = 180
        b_all = 270
        zorder_big = 2
        zorder_small = 1
    if metric == 'gamma':
        ax.set_xlim(-0.3,0.3)
        ax.set_ylim(0,12)
        ymax = 18
        b_opt = 100
        b_true = 200
        b_all = 400
        zorder_big = 2
        zorder_small = 1
    if metric == 'gamma_new':
        ax.set_xlim(-0.3,0.3)
        ax.set_ylim(0,12)
        ymax = 100
        b_opt = 100
        b_true = 200
        b_all = 400
        zorder_big = 2
        zorder_small = 1
    # metric_all = [entry for entry in metric_all if entry<np.mean(metric_all) + 3.5*np.std(metric_all)]
    # metric_true = [entry for entry in metric_true if entry<np.mean(metric_true) + 3.5*np.std(metric_true)]
    metric_optima = [entry for entry in metric_optima if entry<np.mean(metric_optima) + 3.5*np.std(metric_optima)]
    if not is_icsd:
        ax.hist(metric_all, bins = b_all, color = COLORS['orange'], alpha = 0.9, label= 'all', density=True, zorder = 1, edgecolor = "darkorange",)
        ax.hist(metric_true, bins = b_true, color = COLORS['sienna'], alpha = 0.5, label = 'true', density=True, zorder = 3, edgecolor = "saddlebrown",)
        ax.vlines(x = [np.mean(metric_all), np.mean(metric_optima), np.mean(metric_true)], colors = [COLORS['orange'], 
                    COLORS['lightskyblue'], COLORS['sienna']], ymin = 0, ymax = ymax, linestyles='--', linewidth = 2, zorder = 2)
        #ax.vlines(x = [np.mean(metric_optima), np.mean(metric_true)], colors = [ 
                    # COLORS['lightskyblue'], COLORS['sienna']], ymin = 0, ymax = ymax, linestyles='--', linewidth = 2, zorder = 2)
    ax.hist(metric_optima, bins = b_opt, color = COLORS['lightskyblue'], alpha = 0.7, label = 'optimum', density=True, zorder = zorder_big, edgecolor = "darkblue",)
    
    handles1, labels1 = ax.get_legend_handles_labels()
    ax.legend(handles1, labels1, loc='upper right')
    ax.tick_params(left = False)
    ax.set_yticks([])
    # ax2.tick_params(left = False)
    # ax2.legend()
    # plt.legend(handles = [ax1,ax2])
    # return min_data
        #combine all components for one function

def metric_cdf(ax, metric,data):
    data_optimum, data_true, data_all = get_metric_data(metric, data)
    label_colors = {'All': 'orange', 'True':'sienna', 'Optimum':'lightskyblue'}
    for data, label in zip([data_all, data_true, data_optimum], ['All', 'True', 'Optimum']):
    # for data, label in zip([data_true, data_optimum], ['True', 'Optimum']):
        sorted_data = np.sort(data)
        color = label_colors[label]
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        ax.plot(sorted_data, cdf, label=f'{label} CDF', linewidth=2, color = color)

    # Compute Wasserstein Distance
    wd_true_optimum = wasserstein_distance(data_true, data_optimum)
    wd_all_optimum = wasserstein_distance(data_all, data_optimum)
    print(wd_true_optimum, wd_all_optimum)
    ax.set_xlim(-0.5,0.3)
    # Annotate distances
    # ax.set_title(f'CDF Comparison \nW(Optimum, True) = {wd_true_optimum:.4f}, W(Optimum, All) = {wd_all_optimum:.4f}')
    ax.set_xlabel(f'{metric}')
    # ax.set_xlabel('$\Gamma$')
    ax.set_ylabel('Cumulative Probability')
    ax.legend()
    ax.grid()

def metric_hexbin_parity(ax, metric, data):

    # sns.set_theme(style="ticks", rc=set_rc_params())
    m_opt, m_true, m_all = get_metric_data_same_size(metric, data)
    if metric == 'c1':
        gridsize = (90,60)
        vmax = 50
    elif metric == 'c2':
        gridsize = (150,180)
        vmax = 30
    elif metric == 'energy':
        gridsize = (100,33)
        vmax = 30
    elif metric == 'gamma':
        gridsize = (100,80)
        vmax = 50
    elif metric == 'gamma_new':
        gridsize = (100,80)
        vmax = 50
    # Create Hexbin Parity Plot
    hb = ax.hexbin(m_opt, m_true, gridsize=gridsize, cmap='managua', mincnt=1, vmin=0, vmax = vmax)

    # Add Colorbar
    cb = plt.colorbar(hb, ax=ax)
    cb.set_label('Density')


    ax.plot([-10, 10], [-10, 10], linestyle = '--', color = 'red', lw=2, label="Parity Line")
    if metric == 'c1':
        ax.set_xlim(-0.3,0.3)
        ax.set_ylim(-0.3, 0.3)
    elif metric == 'c2':
        ax.set_xlim(-0.01,0.3)
        ax.set_ylim(-0.01,0.3)
    elif metric == 'energy':
        ax.set_xlim(-0.5,0.3)
        ax.set_ylim(-0.5,0.3)
    elif metric == 'gamma':
        ax.set_xlim(-0.3,0.3)
        ax.set_ylim(-0.3,0.3)
    elif metric == 'gamma_new':
        ax.set_xlim(-0.3,0.3)
        ax.set_ylim(-0.3,0.3)
    

    # Labels and Title
    ax.set_xlabel(f'Optimum {metric}')
    ax.set_ylabel(f'Actual {metric}')
    # ax.set_xlabel('Minimum $\Gamma$')
    # ax.set_ylabel('True $\Gamma$')

    # ax.set_title('Hexbin Parity Plot')
    ax.set_aspect('equal')
    # ax.legend()

def all_metric_plots(metric, data):
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    data_list = data_entries(data, solids_data=solids_data)
# Call the plotting functions for each subplot
    metric_histogram(axes[0], metric, data_list)
    metric_cdf(axes[1], metric, data_list)
    metric_hexbin_parity(axes[2], metric, data)

    # Adjust layout
    # fig.suptitle(title)
    plt.tight_layout()
    plt.show()
    return





            
def main():
    solids_data = read_json(os.path.join(DATA_DIR, 'solids_data.json'))
    data_tm = get_formatted_data('tm_closed_metrics.json', optimize = True) 
    # data_tm = get_formatted_data('march_all_metrics.json', optimize = True)
    data_icsd = get_formatted_data('icsd_all_metrics.json', is_icsd=True, optimize = False)
    tm_list = data_entries(data_tm, solids_data=solids_data, filename='tm_data_list.json')
    print('tm_list_done')
    icsd_list = data_entries(data_icsd, solids_data=solids_data, filename='icsd_data_list.json')
    tm_list_opt = data_entries(data_tm, minima=True, solids_data=solids_data)
    icsd_list_opt = data_entries(data_icsd, minima=True, solids_data=solids_data)
    icsd_list_opt = [entry[0] for entry in icsd_list_opt]
    new_materials = read_json(os.path.join(DATA_DIR, '250324_new_materials.json'))
    return solids_data, data_tm, data_icsd, tm_list, icsd_list, tm_list_opt, icsd_list_opt, new_materials

if __name__ == "__main__":
    solids_data, data_tm, data_icsd, tm_list, icsd_list, tm_list_opt, icsd_list_opt, new_materials = main() 
