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
from solidstatesynth.gen.rxn_metrics import PrecursorSet
from scipy.stats import binned_statistic
from mpl_toolkits.axes_grid1 import make_axes_locatable


set_rc_params()
COLORS = CSS4_COLORS # or some other palette
DATA_DIR = "/user/specifies/data/path"


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
    l = []
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
                        for x in d_new:
                            # print(x)
                            x['temperature'] = int(float(temperature))
                            # x['true_rxn'] = False
                    d_new_all[key][str(int(float(temperature)))] = d_new
            else:
                l.append(key)
    print('reformatted',l)
    return d_new_all


def data_filtered(data_r, tm_precursors, ternary = True):
    l = []
    data_new = {key:{str(int(float(t))):[] for t in list(data_r[key].keys())} for key in data_r}
    for key in data_r:
        if ternary:
            if CompTools(key).n_els != 3:
                # print(key)
                continue
        if all([a not in CompTools(key).els for a in ['C', 'N', 'H']]):
            for temperature in data_r[key]:
                for entry in data_r[key][temperature]:
                    rxn_dict = get_reaction_dict_from_string(entry['rxn'])
                    species = rxn_dict['reactants'] + rxn_dict['products']
                    species = [s for s in species if CompTools(s).clean != 'O2']
                    # if all([5>CompTools(s).n_els > 1 for s in species]):
                    if all([7>CompTools(s).n_els for s in species]):

                        if len(rxn_dict['products'])>1:
                            entry['byproducts'] = True
                        else:
                            entry['byproducts'] = False
                        if tm_precursors:
                            if all([CompTools(i).clean in tm_precursors for i in rxn_dict['reactants']]):
                                data_new[key][str(int(float(temperature)))].append(entry)
                            # else:
                            #     print(entry)
                        else:
                            data_new[key][str(int(float(temperature)))].append(entry)

        data_new[key] = {temperature: data_new[key][temperature] 
                            for temperature in data_new[key] if data_new[key][temperature]}
        if not data_new[key]:
            l.append(key)
    # print(l)
    data_new = {key: data_new[key] for key in data_new if data_new[key]}
    return data_new

def data_gamma_new(data_rf, w1, w2, w3):
    for key in data_rf:
        for temp in data_rf[key]:
            for entry in data_rf[key][temp]:
                entry['gamma_new'] = w1*entry['c1'] + w2*entry['c2'] + w3*entry['energy']
    return data_rf


def get_metric_data(metric, data_list, is_icsd = False):
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
        else:
            metric_optima.append(entry[metric])
    return metric_optima, metric_true, metric_all

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


def gamma_comp(weights, entries):
    w1, w2, w3 = weights
    all_diffs = []
    
    for entry in entries:
        true_rxn = [r for r in entry if r['true_rxn']]
        all_true_gamma = [w1 * r['c1'] + w2 * r['c2'] + w3 * r['energy'] for r in true_rxn]
        ### CHRIS EDITS
        if all_true_gamma:
            all_true_gamma = [min(all_true_gamma)]
        else:
            continue
        gamma_opt = min(r['c1'] * w1 + r['c2'] * w2 + r['energy'] * w3 for r in entry)
        diff_list = [(t - gamma_opt)**2 for t in all_true_gamma]
        all_diffs.extend(diff_list)
    # print('med',np.median(all_diffs))
    
    return np.mean(all_diffs)

def optimize_gamma(entries):
    constraints = ({'type': 'eq', 'fun': lambda w: np.sum(w) - 1})  # Constraint: w1 + w2 + w3 = 1
    bounds = [(0.0, 0.99)] * 3  # Ensure w1, w2, w3 are not 0 or 1
    
    # Define the grid of initial values (sum must be 1)
    grid_values = np.linspace(0.01, 0.99, 10)  # Adjust grid density as needed
    initial_guesses = [p for p in product(grid_values, repeat=3) if np.isclose(sum(p), 1.0, atol=1e-2)]
    # print('i',initial_guesses)
    best_result = 1e6
    best_w = None

    for guess in initial_guesses:
        result = minimize(gamma_comp, guess, args=(entries,), bounds=bounds, constraints=constraints)
        # print(guess, result.x, result.fun)
        if result.fun < best_result:
            best_result = result.fun
            best_w = result.x
            # print(best_w, best_result)

    return best_w  # Optimized w1, w2, w3


def get_formatted_data(file_name, tm_precursors, tm_rxns = None, optimize = True, is_icsd = False, ternary = True):
    # count = 0
    data = read_json(os.path.join(DATA_DIR_metrics, file_name))
    if type(data) == list:
        data = {key: entry[key] for entry in data for key in entry}
    data = data_reformatted(data, is_icsd)
    # if 'Cu1Nb2O6' in data:
        # print('yes ref')
    data = data_filtered(data, tm_precursors=tm_precursors, ternary=ternary)
    # if 'Cu1Nb2O6' in data:
    #     print('yes filter')
    if tm_rxns:
        data_i = find_true_rxns(tm_rxns, data)
        print(len(data_i))
    else:
        data_i = data
    # if 'Cu1Nb2O6' in data_i:
    #     print(data_i['Cu1Nb2O6'])
    if is_icsd:
        return data_gamma_new(data, 0.562,0.280,0.158)
    data_n = {}
    for key in data_i:
        for temperature in data_i[key]:
            if any([entry['true_rxn'] == True for entry in data_i[key][temperature]]):
                if key not in data_n:
                    data_n[key] = {}
                data_n[key][temperature] = data_i[key][temperature]
                # count += 1
    # print(count)
    if optimize:
        data_list = []
        for key in data_n:
            # print(key)
            for temp in data_n[key]:
                if data_n[key][temp]:
                    data_list.append(data_n[key][temp])
        # print(data_list[0:2])
        w1, w2, w3 = optimize_gamma(data_list)
        print('opt', w1,w2,w3)
    else:
        w1, w2, w3 =0.562,0.280,0.158
    # print(w1,w2,w3)
    data_new = data_gamma_new(data_n, w1,w2,w3)
    # print(count)
    return data_new

def get_target_stability(target, temperature, stability):
    if target not in stability:
        # print(target)
        # print(target)
        return None
    if temperature not in stability[target]:
        if temperature == '1073':
            if '1073.15' in stability[target]:
                return stability[target]['1073.15']
            else:
                return None
        if temperature == '1073.15':
            if '1073' in stability[target]:
                return stability[target]['1073']
            else:
                return None
        return None
    return stability[target][temperature]

def data_entries(data, stability, minima = False, gamma = 'gamma_new'):
    # p = os.path.join(DATA_DIR,filename)
    # if os.path.exists(p):
    #     return read_json(p)['data']
    data_entry_list = []
    for key in data:
        if minima == False:
            for temperature in data[key]:
                min_gamma = min([entry[gamma] for entry in data[key][temperature]])
                min_entry = [entry for entry in data[key][temperature] if entry[gamma] == min_gamma][0]
                min_entry['opt_rxn'] = True
                for entry in data[key][temperature]:
                    if temperature == 'null':
                        temperature = 1073.15
                    entry['temperature'] = int(float(temperature))
                    entry['target'] = CompTools(key).clean
                    entry['E_d'] = get_target_stability(CompTools(key).clean, temperature, stability)
                    if 'opt_rxn' not in entry:
                        entry['opt_rxn'] = False
                    data_entry_list.append(entry)

        else:
            mins = []
            for temperature in data[key]:
                min_gamma = [entry[gamma] for entry in data[key][temperature]]
                if not min_gamma:
                    print(key, temperature)
                    continue
                else:
                    min_gamma = min(min_gamma)
                min_entry = [entry for entry in data[key][temperature] if entry[gamma] == min_gamma][0]
                if temperature == 'null':
                    temperature = 1073.15
                min_entry['temperature'] = int(float(temperature))
                min_entry['target'] = CompTools(key).clean
                min_entry['E_d'] = get_target_stability(CompTools(key).clean, temperature, stability)

                mins.append(min_entry)
            if mins:
                min_min = min([entry[gamma] for entry in mins])
                min_entry = [entry for entry in mins if entry[gamma] == min_min]
                data_entry_list.extend(min_entry)
    # f = write_json({'data':data_entry_list},p)
    return data_entry_list



# Example histogram data (replace with actual histograms)

def metric_histogram(ax, metric, data, icsd = None):
    metric_optima, metric_true, metric_all = get_metric_data(metric, data)
    if icsd:
        icsd_metric_optima = get_metric_data(metric, icsd, is_icsd=True)[0]
        # print(np.mean(icsd_metric_optima))
    # fig, ax1 = plt.subplots()
    # ax.set_xlabel(f"{metric}")
    ax.set_xlabel("$\Gamma$ (eV/atom)")
    if metric == 'c1':
        ax.set_xlim(-1,1)
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
        b_opt = 200
        b_true = 200
        b_all = 400
        zorder_big = 2
        zorder_small = 1
    if metric == 'gamma_new':
        ax.set_xlim(-0.3,0.4)
        ax.set_ylim(0,12)
        ymax = 100
        # b_opt = 150
        b_opt = 200
        b_true = 300
        b_all = 225
        zorder_big = 2
        zorder_small = 1

    if not icsd:
        # ax.hist(metric_all, bins = b_all, color = COLORS['orange'], alpha = 0.9, label= 'all', density=True, zorder = 1, edgecolor = "darkorange",)
        ax.hist(metric_true, bins = b_true, color = COLORS['sienna'], alpha = 0.5, label = 'observed', density=True, zorder = 3, edgecolor = "saddlebrown",)
        ax.hist(metric_optima, bins = b_opt, color = COLORS['lightskyblue'], alpha = 0.7, label = 'optimum', density=True, zorder = zorder_big, edgecolor = "darkblue",)
        ax.vlines(x=0.088, ymin = -0.1, ymax = 20, color = 'maroon', linestyle = '--', linewidth = 2, zorder = 10)

        # ax.vlines(x = [np.mean(metric_true), np.mean(metric_optima)], colors = ['red', 'blue'], ymin = 0, ymax = 15, linestyles='--', linewidth = 2, zorder = 12)

    else:
        ax.hist(metric_true, bins = b_true, color = COLORS['sienna'], alpha = 0.5, label = 'TM true rxns', density=True, zorder = 3, edgecolor = "saddlebrown",)
        ax.hist(metric_optima, bins = b_opt, color = COLORS['lightskyblue'], alpha = 0.7, label = 'TM optimum rxns', density=True, zorder = zorder_big, edgecolor = "darkblue",)
        ax.hist(icsd_metric_optima, bins = 200, color = COLORS['forestgreen'], alpha = 0.7, label = 'ICSD', density=True, zorder = 10, edgecolor = "darkgreen",)
        print(np.mean(icsd_metric_optima))
        # ax.vlines(x = [np.mean(metric_true), np.mean(metric_optima)], colors = ['blue', 'green'], ymin = 0, ymax = 15, linestyles='--', linewidth = 2, zorder = 12)
    handles1, labels1 = ax.get_legend_handles_labels()
    ax.legend(handles1, labels1, loc='upper left')
    ax.tick_params(left = False)
    ax.set_xticks([-0.2,0,0.2])
    ax.set_yticks([])
    # ax.set_ylabel('Counts')

def metric_cdf(ax, metric,data, icsd = None):
    data_optimum, data_true, data_all = get_metric_data(metric, data)

    label_colors = {'all': 'orange', 'observed':'sienna', 'optimum':'lightskyblue', 'ICSD':'forestgreen'}
    if not icsd:
        for data, label in zip([data_true, data_optimum], ['observed', 'optimum']):
        # for data, label in zip([data_true, data_optimum], ['True', 'Optimum']):
            sorted_data = np.sort(data)
            color = label_colors[label]
            cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
            ax.plot(sorted_data, cdf, label=f'{label} CDF', linewidth=2, color = color)
    else:
        label_colors = {'all': 'orange', 'observed':'sienna', 'optimum':'lightskyblue', 'ICSD optimum rxns':'forestgreen'}

        icsd_data = get_metric_data(metric,icsd, is_icsd = True)[0]
        for data, label in zip([icsd_data, data_true, data_optimum], ['ICSD optimum rxns', 'TM true rxns', 'TM optimum rxns']):
        # for data, label in zip([data_true, data_optimum], ['True', 'Optimum']):
            sorted_data = np.sort(data)
            color = label_colors[label]
            cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
            ax.plot(sorted_data, cdf, label=f'{label} CDF', linewidth=2, color = color)

    # Compute Wasserstein Distance
    wd_true_optimum = wasserstein_distance(data_true, data_optimum)
    wd_all_optimum = wasserstein_distance(data_all, data_optimum)

    if icsd:
        wd_icsd_optimum = wasserstein_distance(icsd_data, data_optimum)
        print('i,opt',wd_icsd_optimum)
    print('t,opt', wd_true_optimum)
    print(wd_all_optimum)

    ax.set_xlim(-0.2,0.36)
    # ax.hlines(y=0.9, xmin = -0.5, xmax = 0.5, color = 'black', linestyle = '--', linewidth = 2)
    # ax.hlines(y=0.68, xmin = -0.5, xmax = 0.5, color = 'blue', linestyle = '--')
    # ax.hlines(x=0.2, ymin = -0.1, ymax = 1, color = 'red', linestyle = '--')
    ax.vlines(x=0.088, ymin = -0.1, ymax = 1, color = 'maroon', linestyle = '--', linewidth = 2)
    # ax.vlines(x=[0.06,0.325], ymin = -0.1, ymax = 1, color = 'red', linestyle = '--')
    ax.set_ylim(-0.01,1.01)

    # Annotate distances
    # ax.set_title(f'CDF Comparison \nW(Optimum, True) = {wd_true_optimum:.4f}, W(Optimum, All) = {wd_all_optimum:.4f}')
    # ax.set_xlabel(f'{metric}')
    ax.set_xlabel('$\Gamma$ (eV/atom)')
    ax.set_ylabel('Cumulative Probability')
    ax.legend(loc = 'lower right')
    ax.grid()
    ax.set_yticks([0,0.5,1.0])

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
    # cb = plt.colorbar(hb, ax=ax, ticks = [0,50],shrink = 0.75)
    # cb.set_label('Number of targets')


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
    ax.set_xlabel(f'Opt. $\Gamma$ (eV/atom)')
    ax.set_ylabel(f'Obs. $\Gamma$ (eV/atom)')
    # ax.set_xlabel('Minimum $\Gamma$')
    # ax.set_ylabel('True $\Gamma$')

    # ax.set_title('Hexbin Parity Plot')
    # ax.set_aspect('equal')
    # ax.legend()
    ax.set_yticks([-0.2,0,0.2])
    ax.set_xticks([-0.2,0,0.2])
    ax.grid()




def plot_binned_average(ax, data, title='', xlabel='E$_{hull}$ (eV/atom)', ylabel='$\Gamma$ (eV/atom)', 
                        xlim=(0, 0.1), ylim=(-0.2, 0.4), bins=15, score_line=None):
    """
    Create a binned average plot with error bars (standard deviation) and connecting lines.
    """
    # Compute binned means and standard deviations
    data = [entry for entry in data if entry['true_rxn']]
    data = [entry for entry in data if entry['E_d'] != None]
    y = [entry['gamma_new'] for entry in data]
    x = []
    for entry in data:
        if entry['E_d']<0:
            x.append(0)
        else:
            x.append(entry['E_d'])

    bs_avg, edges, _ = binned_statistic(x, y, statistic='mean', bins=bins, range=xlim)
    bs_std, _, _ = binned_statistic(x, y, statistic='std', bins=bins, range=xlim)

    bin_centers = (edges[1:] + edges[:-1]) / 2

    # Plot the binned average with connecting lines and error bars
    # ax.figure(figsize=(6, 4))
    ax.errorbar(bin_centers, bs_avg, yerr=bs_std, fmt='-o', color='blue', ecolor='deepskyblue', 
                 elinewidth=2, capsize=3, label='Binned Average')

    # Optional score line
    if score_line:
        ax.plot(np.linspace(xlim[0], xlim[1], 100), np.ones(100) * score_line,
                 color='red', linestyle='--', linewidth=2, label=f'Score line: {score_line}')

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.legend()
    # ax.set_aspect('equal')
    ax.set_xticks([0,0.05,0.1])
    ax.set_yticks([-0.2,0,0.2, 0.4])
    ax.grid()


def all_metric_plots(metric, data, stability, icsd = None, sz = 30):
    data_list = data_entries(data, stability)

    fig, axes = plt.subplots(1, 4, figsize=(30, 6), constrained_layout = True)
    hb = metric_hexbin_parity(axes[1], metric, data)
    # fig.colorbar(c, ax=axs[1, 1], fraction=0.046, pad=0.04)

    # else:
    #     fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    plot_binned_average(axes[0], data_list)

    metric_histogram(axes[2], metric, data_list, icsd = icsd)
    metric_cdf(axes[3], metric, data_list, icsd=icsd)
    # Adjust layout
    # fig.suptitle(title)
    # plt.tight_layout()
    m_opt, m_true, m_all = get_metric_data_same_size(metric, data)
    # hb = ax.hexbin(m_opt, m_true, gridsize=gridsize, cmap='managua', mincnt=1, vmin=0, vmax = vmax)
    fig.colorbar(mappable= hb, ax=axes[1], fraction=0.046, pad=0.04, cmap = 'managua')
    plt.rcParams.update({
        'font.size': 0.9*sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.9*sz
    })    
    # plt.subplots_adjust(right = 0.4)

    plt.show()
    return

def tm_data():
    tm_chemsys_data = read_json(os.path.join(DATA_DIR, 'new_tm_chemsys.json'))
    tm = {}
    count = 0 
    for key in tm_chemsys_data:
        d = tm_chemsys_data[key]
        for entry in d:
            if entry['atmosphere'] not in ['air', None]:
                continue
            if 'nan' in entry['rxn']:
                continue
            if '-' in entry['rxn']:
                continue
            rxn_dict = get_reaction_dict_from_string(entry['rxn'])
            try:
                precursors = set([CompTools(i).clean for i in rxn_dict['reactants']])
                if 'O1' in precursors:
                    'O1'.replace('O2')
                entry['precursors'] = precursors
            except:
                # print('tm unclean')
                continue
            try:
                products = set([CompTools(i).clean for i in rxn_dict['products']])
            except:
                # print('tm unclean')
                continue
            targ = CompTools(entry['target']).clean
            if len(CompTools(targ).els) == 3:
                if 'O' in targ:
                    if targ not in tm:
                        tm[targ] = []
                        tm[targ].append(entry)
                        count += 1
                    else:
                        td = []
                        if entry['temperature'] == None:
                            e_temp = 1073
                        else:
                            e_temp = entry['temperature']
                        for r in tm[targ]:
                            if not r['temperature']:
                                temp = 1073
                            else:
                                temp = r['temperature']
                            if abs(temp - e_temp)<50:
                                td.append(r)
                        if not td:
                            tm[targ].append(entry)
                            count += 1
                        else:
                            same_p = [t for t in td if t['precursors']==precursors]
                            if not same_p:
                                tm[targ].append(entry)
                                count += 1
            
    print(count)
    return tm

def is_true_rxn(tm_rxn, rxn_metrics_entry):
    rxn_dict = get_reaction_dict_from_string(tm_rxn['rxn'])
    try:
        precursors = set([CompTools(i).clean for i in rxn_dict['reactants']])
    except:
        # print('tm unclean')
        return False
    try:
        products = set([CompTools(i).clean for i in rxn_dict['products']])
    except:
        # print('tm unclean')
        return False
    # print(precursors, products)
    entry_rxn = rxn_metrics_entry['rxn']
    entry_rxn_dict = get_reaction_dict_from_string(entry_rxn)
    try:
        entry_precursors = set([CompTools(i).clean for i in entry_rxn_dict['reactants']])
    except:
        # print('unclean')

        return False
    try:
        entry_products = set([CompTools(i).clean for i in entry_rxn_dict['products']])
    except:
        # print('unclean')

        return False
    if entry_precursors == precursors:
        if entry_products == products:
            return True
    return False

def find_true_rxns(tm_rxns, data_tm):
    af = []
    nt = []
    nk = []
    for key in tm_rxns:
        if key not in data_tm:
            # print('not in data',key)
            nk.append(key)
            continue
        calc_rxns = data_tm[key]
        for entry in tm_rxns[key]:
            entry['matched'] = False
            temp = entry['temperature']
            if temp in ['None',None]:
                temp = 1073
            temp = str(int(float(temp)))
            if temp in calc_rxns:
                for rxn in calc_rxns[temp]:
                    if is_true_rxn(entry, rxn):
                        rxn['true_rxn'] = True
                        entry['matched'] = True
                    else:
                        rxn['true_rxn'] = False
                if all([rxn['true_rxn'] == False for rxn in calc_rxns[temp]]):
                    # print('all false', key)
                    af.append(key)
            else:
                ts = [t for t in calc_rxns if abs(float(t)-float(temp))<50]
                # print(ts)
                if ts:
                    for rxn in calc_rxns[ts[0]]:
                        if is_true_rxn(entry, rxn):
                            rxn['true_rxn'] = True
                            entry['matched'] = True
                        else:
                            rxn['true_rxn'] = False
                    if all([rxn['true_rxn'] == False for rxn in calc_rxns[ts[0]]]):
                        # print('all false', key)
                        af.append(key)
                else:
                    nt.append(key)
                    # print(ts)
                    # print(calc_rxns.keys(), key, temp)
    return data_tm

def get_id_from_formula(formula, solids_data):
    return solids_data[formula]['material_id']

def synth_scores(synth, solids_data, data_tm):
    scores = {}
    formulas = list(set([f for f in data_tm]))
    for f in formulas:
        opt_gammas = [entry['gamma_new'] for temp in data_tm[f] for entry in data_tm[f][temp] if entry['opt_rxn']]
        if not opt_gammas:
            continue
        else:
            opt_gamma = min(opt_gammas)
        if f not in solids_data:
            print(f)
            continue
        mpid = solids_data[f]['material_id']
        Ed = solids_data[f]['energy_above_hull']
        if mpid not in synth:
            continue
        score = synth[mpid]
        scores[f] = {'Ed':Ed, 'pred':score, 'opt_gamma':opt_gamma}
    return scores

            
def main():
    solids_data = read_json(os.path.join(DATA_DIR, 'solids_data.json'))
    sd = read_json(os.path.join(DATA_DIR, '241119_mp_gd.json'))['data']
    sd = {CompTools(entry['formula']).clean: entry for entry in sd}
    solids_data = {CompTools(k).clean: solids_data[k] for k in solids_data}
    tm_precursors = read_json(os.path.join(DATA_DIR, 'tm_precursors.json'))['data']
    tm_rxns = tm_data()
    stability = read_json(os.path.join(DATA_DIR, 'stability_data_tm.json'))
    synth_pred = read_json(os.path.join(DATA_DIR_gen, 'mpid_predictions.json'))
    stability_i = read_json(os.path.join(DATA_DIR, 'stability_data_tm_new_i.json'))
    stability = {CompTools(key).clean:stability[key] for key in stability}
    stability_i = {CompTools(key).clean:stability_i[key] for key in stability_i}
    stability.update(stability_i)
    len(stability)
    s_new = {}
    for key in stability:
        s_new[key] = {}
        for temp in stability[key]:
            if temp == 'null':
                t_new = '1073'
            else:
                t_new = str(int(float(temp)))
            s_new[key][t_new] = stability[key][temp]
    stability = s_new
    data_tm = get_formatted_data('250506_tm_closed_metrics.json', tm_rxns=tm_rxns,tm_precursors=None, optimize = True, ternary=True)
    # print(len(data_tm))
    data_tm = {key:data_tm[key] for key in data_tm if CompTools(key).n_els==3}
    data_gen = get_formatted_data('gen_closed_metrics.json',tm_precursors=None, is_icsd=True)
    gen_list = data_entries(data_gen, stability=stability, minima=True)
    print(len(data_tm))
    print('tm')
    tm_list = data_entries(data_tm, stability=stability)
    tm_true = [entry for entry in tm_list if entry['true_rxn']]
    tm_opt = [entry for entry in tm_list if entry['opt_rxn']]
    print(len(tm_true), 'true')
    tm_list_opt = data_entries(data_tm, stability=stability, minima=True)
    pred_dict = synth_scores(synth_pred,sd, data_tm)
    # write_json(pred_dict, os.path.join(DATA_DIR,'tm_predictions.json'))
    return solids_data, stability, stability_i, synth_pred, data_tm, tm_list, tm_list_opt, tm_true, tm_opt, pred_dict, data_gen, gen_list

if __name__ == "__main__":
    solids_data, stability, stability_i, synth_pred, data_tm, tm_list, tm_list_opt, tm_true, tm_opt, pred_dict, data_gen, gen_list = main() 
