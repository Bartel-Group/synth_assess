import numpy as np
from scipy.stats import entropy

def jensen_shannon_divergence(p, q):
    """
    Compute the Jensen-Shannon Divergence between two probability distributions.
    """
    p = np.array(p) + 1e-10  # Avoid zeros
    q = np.array(q) + 1e-10  # Avoid zeros
    m = 0.5 * (p + q)  # Mixture distribution

    return 0.5 * (entropy(p, m) + entropy(q, m))  # JS Divergence

# Compute histograms and normalize

def metric_histogram(metric, data):
    metric_optima = []
    metric_true = []
    metric_all = []
    # min_data = []
    for key in data:
        for temperature in data[key]:
            temp_data = data[key][temperature]
            if temp_data:
                all_temp = [entry[metric] for entry in temp_data]
                metric_all.extend(all_temp)
                metric_optima.append(min(all_temp))
                # min_data_o = [entry for entry in temp_data if entry[metric] == min(metric_all)]
                metric_true.extend([entry[metric] for entry in temp_data if entry['true_rxn']])
    bins = np.linspace(-0.5, 0.5, 30)
    hist_all, _ = np.histogram(metric_all, bins=bins, density=True)
    hist_true, _ = np.histogram(metric_true, bins=bins, density=True)
    hist_optimum, _ = np.histogram(metric_optima, bins=bins, density=True)


    # Compute JS divergence
    js_true_optimum = jensen_shannon_divergence(hist_true, hist_optimum)
    js_all_optimum = jensen_shannon_divergence(hist_all, hist_optimum)

    print(f"JS Divergence (True vs. Optimum): {js_true_optimum:.4f}")
    print(f"JS Divergence (All vs. Optimum): {js_all_optimum:.4f}")

    if js_true_optimum < js_all_optimum:
        print("The blue histogram is statistically closer to the green histogram than the orange histogram.")
    else:
        print("The orange histogram is closer to the green histogram.")

    # print('true vs. opt', jensen_shannon_divergence(metric_true, metric_optima))
    # print('all vs. opt', jensen_shannon_divergence(metric_all, metric_optima))
    return metric_all, metric_optima, metric_true

import numpy as np
from scipy.stats import entropy

def jensen_shannon_divergence(p, q):
    """Compute the JS divergence between two probability distributions."""
    p = np.array(p) + 1e-10  # Avoid zeros
    q = np.array(q) + 1e-10
    m = 0.5 * (p + q)
    return 0.5 * (entropy(p, m) + entropy(q, m))

def bootstrap_js(all, true, opt, bins=30, n_bootstrap=1000):
    """
    Perform bootstrapping to compute confidence intervals for JS divergence.
    """
    js_diffs = []
    data1 = all
    data2 = true
    data3 = opt
    
    for _ in range(n_bootstrap):
        # Resample data with replacement
        sample1 = np.random.choice(data1, size=len(data1), replace=True)
        sample2 = np.random.choice(data2, size=len(data2), replace=True)
        sample3 = np.random.choice(data3, size=len(data3), replace=True)
        
        # Compute histograms and normalize
        hist1, _ = np.histogram(sample1, bins=bins)
        hist2, _ = np.histogram(sample2, bins=bins)
        hist3, _ = np.histogram(sample3, bins=bins)
        
        hist1 = hist1 / np.sum(hist1)
        hist2 = hist2 / np.sum(hist2)
        hist3 = hist3 / np.sum(hist3)
        
        # Compute JS divergence
        js_true_optimum = jensen_shannon_divergence(hist2, hist3)  # True vs Optimum
        js_all_optimum = jensen_shannon_divergence(hist1, hist3)  # All vs Optimum
        
        # Store the difference
        js_diffs.append(js_all_optimum - js_true_optimum)

    # Compute confidence interval
    lower = np.percentile(js_diffs, 2.5)
    upper = np.percentile(js_diffs, 97.5)
    mean_diff = np.mean(js_diffs)

    print(f"Bootstrapped JS difference (All - True): {mean_diff:.4f}")
    print(f"95% Confidence Interval: [{lower:.4f}, {upper:.4f}]")

    # Interpretation
    if lower > 0:
        print("The difference is statistically significant: 'All' is farther from 'Optimum' than 'True'.")
    elif upper < 0:
        print("The difference is significant in the opposite direction: 'True' is farther than 'All'.")
    else:
        print("The difference is not statistically significant.")
    
    return mean_diff, lower, upper


def view_metric(data):
    energy_optima = []
    c1_optima = []
    gamma_optima = []
    # min_data = []
    for key in data:
        for temperature in data[key]:
            temp_data = data[key][temperature]
            if temp_data:
                min_gamma = min(temp_data, key=lambda d: d["gamma"])
                energy_optima.append(min_gamma['energy'])
                c1_optima.append(min_gamma['c1'])
                gamma_optima.append(min_gamma['gamma'])
                
    # all_stdev = np.std(metric_all)
    fig, ax2 = plt.subplots()
    # ax1.set_xlabel(f"{metric}")
    # ax1.set_xlim(np.mean(metric_all) - 2.5*all_stdev, np.mean(metric_all) + 2*all_stdev)
    # ax1.scatter(gamma_optima, gamma_optima, color = COLORS['orange'], alpha = 0.8, label= 'all')
    # ax1.legend()

    ax2.scatter(c1_optima, gamma_optima, color = COLORS['blue'], alpha = 0.8, label = 'c1')
    ax2.scatter(energy_optima, gamma_optima, color = COLORS['green'], alpha = 0.8, label = 'energy')
    # ax2.plot([-0.5,0.5],[-0.5,0.5], linestyle = '--')
    # ax2.set_ylim([0,260])
    # handles1, labels1 = ax1.get_legend_handles_labels()
    # handles2, labels2 = ax2.get_legend_handles_labels()
    # handles = handles1 + handles2
    # labels = labels1 + labels2
    ax2.legend()
    ax2.set_xlim(-0.5,0.5)
    ax2.set_ylim(-0.5,0.5)
    ax2.set_xlabel('optimum of metric')
    ax2.set_ylabel('optimum gamma')
    # # Add a single legend
    # ax1.legend(handles, labels, loc='upper right')
    # ax1.tick_params(left = False)
    # ax2.tick_params(left = False)
    # ax2.legend()
    # plt.legend(handles = [ax1,ax2])
    plt.show()
    # return min_data
        #combine all components for one function

def compare_metrics(data):
    metric_optima_c1 = []
    metric_true_c1 = []
    metric_all_c1 = []
    metric_optima_c2 = []
    metric_true_c2 = []
    metric_all_c2 = []
    metric_optima_gamma = []
    metric_true_gamma = []
    metric_true_energy = []
    min_c1s = []
    min_c2s = []
    min_energy = []
    metric_optima_energy = []
    # min_data = []
    for key in data:
        for temperature in data[key]:
            temp_data = data[key][temperature]
            if temp_data:
                min_gamma = min(temp_data, key=lambda d: d["gamma"])
                all_temp_c1 = [entry['c1'] for entry in temp_data]
                min_c1 = min(all_temp_c1)
                all_temp_c2 = [entry['c2'] for entry in temp_data]
                energy_all = [entry['energy'] for entry in temp_data]
                metric_all_c1.extend(all_temp_c1)
                metric_all_c2.extend(all_temp_c2)
                metric_optima_c1.append(min_gamma['c1'])
                min_c1s.append(min_c1)
                min_c2s.append(min(all_temp_c2))
                min_energy.append(min(energy_all))
                metric_optima_c2.append(min_gamma['c2'])
                metric_optima_gamma.append(min_gamma['gamma'])
                metric_optima_energy.append(min_gamma['energy'])
                # min_data_o = [entry for entry in temp_data if entry[metric] == min(metric_all)]
                metric_true_c1.extend([entry['c1'] for entry in temp_data if entry['true_rxn']])
                metric_true_c2.extend([entry['c2'] for entry in temp_data if entry['true_rxn']])
                metric_true_gamma.extend([entry['gamma'] for entry in temp_data if entry['true_rxn']])
                metric_true_energy.extend([entry['energy'] for entry in temp_data if entry['true_rxn']])
                # min_data.extend(min_data_o)
                # print(len(min_data))
    # norm_factor_true = int(len(metric_all)/3*len(metric_true))
    # norm_factor_opt = int(len(metric_all)/2*len(metric_optima))
    # plt.scatter(metric_all_c1, metric_all_c2, color = COLORS['orange'], label = 'all')
    # plt.hist(metric_optima_gamma, label = 'optimum energy')

    c1_diff = [abs(min_c2s[i] - metric_optima_c2[i]) for i in range(len(metric_optima_c1))]
    # print(len([entry for entry in c_diff if entry == 0]))
    # print(np.mean([entry for entry in c_diff if entry != 0]))
    c2_diff = [abs(min_c1s[i] - metric_optima_c1[i]) for i in range(len(metric_optima_c1))]
    c1_opt = [i for i in range(len(c1_diff)) if c1_diff[i] == 0]
    c2_opt = [i for i in range(len(c2_diff)) if c2_diff[i] == 0]
    c2_d = [c2_diff[i] for i in range(len(c1_diff)) if c1_diff[i] == 0 if c2_diff[i] != 0] 
    c2_d = [entry for entry in c2_d if entry < 0.3]
    print(np.mean(c2_d))
    print(c2_d)
    print(len(c2_d))
    c1_d = [c1_diff[i] for i in range(len(c2_diff)) if c2_diff[i] == 0 if c1_diff[i] != 0]
    c1_d = [entry for entry in c1_d if entry < 0.3]

    print(len(c1_d))
    print(np.mean(c1_d))
    print(len(c1_opt))
    print(len(c2_opt))

    # print(len(both_opt))
    # print(len([entry for entry in c_diff if entry == 0]))
    # print(np.mean([entry for entry in c_diff if entry != 0]))
    # print(len(c_diff))
    plt.scatter(c1_diff, c2_diff)
    # plt.scatter(metric_optima_gamma, metric_optima_energy)
                # c = metric_optima_c2, cmap = 'viridis', vmin = 0, vmax = 0.5, label = 'c1')
    # plt.colorbar()
    # for i, item in enumerate(metric_optima_gamma):
    #     if metric_
    # a, b = np.polyfit(metric_optima_c1, metric_optima_gamma, 1)
    # y_hyp = [entry*a + b for entry in metric_optima_c1]
    # plt.plot(metric_optima_c1, y_hyp, linestyle = '--')
    # plt.scatter(metric_optima_gamma, metric_optima_c2, label = 'c2')
    plt.ylabel('c2 diff')
    plt.xlabel('c1 diff')
    plt.xlim(-0.01, 0.2)
    plt.ylim(-0.01 ,0.2)
    # plt.ylim(-1, 0.5)
    # plt.plot([0, 1], [0,1], linestyle = '--')
    # plt.xlim(-0.01,0.5)
    # plt.ylim(0, 0.5)
    # plt.ylim(0, 1)
    plt.legend()
    plt.show()
    # return c_diff
        #combine all components for one function"
