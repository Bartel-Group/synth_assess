import numpy as np
from matplotlib.patches import PathPatch
from matplotlib.path import Path
import matplotlib.pyplot as plt
import matplotlib as mpl


# Helper functions for figure 1

def plot_hull(ax, target, ylabel = '$\Delta$$\it{G}$$_{rxn}$ ' + '(eV/atom)', p_sz = 300):
    """ 
    helper function for Figure 1 plotting
    """
   
    if target == 'LaNb$_3$O$_9$':
        precs = ['Nb$_2$O$_5$','La$_2$O$_3$']
        coords = [[0.125, 0.5, 0.75],[-0.051,-0.204,-0.148]]
        target_coords = [0.25,-0.093]
        color2 = 'red'
        color1 = 'red'
        ax.text(target_coords[0],target_coords[1]+0.01,target, color = 'red')
    elif target == 'LaGaO$_3$':
        precs = ['Ga$_2$O$_3$','La$_2$O$_3$']
        coords = [[0.375,0.667],[-0.081, -0.1]]
        target_coords = [0.5,-0.043]
        color2 = 'royalblue'
        color1 = 'cornflowerblue'
        ax.text(target_coords[0]+0.02,target_coords[1]+0.01,target, color = 'royalblue')
    for i in range(len(coords[0])):
        if i == 0:
            ax.plot((0,coords[0][0]), (0, coords[1][0]), color = 'black', linewidth = lw)
        if i == len(coords[0])-1:
            ax.plot((coords[0][i],1), (coords[1][i],0), color = 'black', linewidth = lw)
        else:
            ax.plot((coords[0][i],coords[0][i+1]), (coords[1][i],coords[1][i+1]), color = 'black', linewidth = 2.5)
    ax.scatter([0]+coords[0]+[1],[0]+coords[1]+[0], s = p_sz, edgecolor = color1, color = 'black', zorder = 4)
    ax.scatter([target_coords[0]],[target_coords[1]], s = p_sz, color = color2, edgecolor='black', zorder = 4)
    # ax.grid(True)
    ax.text(0.02,0.01, precs[0])
    ax.text(0.71,0.01, precs[1])
    ax.hlines(y = 0, xmin = 0, xmax = 1, color = 'black', linewidth = 1.5)
    ax.set_xlabel('$\it{x}$ in (1 - $\it{x})$'+precs[0] + ' + ' + '$\it{x}$' + precs[1])
    ax.set_ylabel(ylabel)
    ax.set_ylim(-0.22,0.05)
    ax.set_xlim(0,1)
    ax.set_xticks([0,0.25,0.5,0.75,1])
    ax.set_xticklabels(['0','','0.5','','1.0'])

def plot_rolling_stats(ax, Ehull, Gamma, window=0.008, step=0.0025, 
                       color1='skyblue', scatter_alpha=0.3,
                       xlim=(0, 0.1), ylim=(-0.16, 0.32)):
    # Convert inputs to numpy arrays in case they're lists 
    """
    helper function for Figure 1 plotting

    """
    Ehull = np.array(Ehull)
    Gamma = np.array(Gamma)

    # Filter to desired plotting range
    mask = (Ehull >= xlim[0]) & (Ehull <= xlim[1]) & (Gamma >= ylim[0]) & (Gamma <= ylim[1])
    Ehull_filtered = Ehull[mask]
    Gamma_filtered = Gamma[mask]

    # Prepare DataFrame for rolling stats
    df = pd.DataFrame({'Ehull': Ehull_filtered, 'Gamma': Gamma_filtered}).sort_values('Ehull')
    
    x_vals = np.arange(xlim[0], xlim[1], step)
    
    medians, p25s, p75s, x_centers = [], [], [], []


    for x in x_vals:
        window_mask = (df['Ehull'] >= x) & (df['Ehull'] < x + window)
        if window_mask.sum() > 0:
            subset = df[window_mask]['Gamma']
            medians.append(subset.median())
            p25s.append(subset.quantile(0.25))
            p75s.append(subset.quantile(0.75))
            x_centers.append(x + window / 2)

    x_centers = np.array(x_centers)
    medians = np.array(medians)
    p25s = np.array(p25s)
    p75s = np.array(p75s)

    # Extend rolling stats to full xlim
    if len(x_centers) > 0:
        # Left edge
        if x_centers[0] > xlim[0]:
            x_centers = np.insert(x_centers, 0, xlim[0])
            medians   = np.insert(medians,   0, medians[0])
            p25s      = np.insert(p25s,      0, p25s[0])
            p75s      = np.insert(p75s,      0, p75s[0])

        # Right edge
        if x_centers[-1] < xlim[1]:
            x_centers = np.append(x_centers, xlim[1])
            medians   = np.append(medians,   medians[-1])
            p25s      = np.append(p25s,      p25s[-1])
            p75s      = np.append(p75s,      p75s[-1])


    # Scatter plot of raw data
    ax.scatter(Ehull_filtered, Gamma_filtered, alpha=scatter_alpha, color='gray', label='Individual recipe', s = 100)

    # Shaded area: 25th–75th percentile
    ax.fill_between(x_centers, p25s, p75s, color=color1, alpha=0.5, label='25th–75th percentile')

    # Dashed median line
    ax.plot(x_centers, medians, linestyle='--', color='black', label='Rolling median', linewidth = 2.5)

    # Set axis limits and labels
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_xticklabels(['0.0','0.02','0.04','0.06','0.08','0.1'])
    ax.set_xticks([0,0.02,0.04,0.06, 0.08, 0.1])
    ax.set_yticks([-0.1,0,0.1,0.2,0.3])
    ax.set_xlabel("$\it{E}$$_{hull}$ (eV/atom)")
    ax.set_ylabel("$\Gamma_{obs}$ (eV/atom)")
    ax.legend(frameon = True, framealpha = 1, loc = 'best')
    # ax.text(x=0, y = -0.16, s = sz, text = '0.0')     
    # ax.grid(True)

# Figure 2

def metric_cdf(ax, tm_entries):
    """
    helper function for plotting figure 2
    """
    metric_optima = []
    metric_true = []
    data_optimum = []
    data_true = []
    for entry in tm_entries:
        if entry['true_rxn']:
            metric_true.append(entry['gamma_new'])
            data_true.append(entry)
        if entry['opt_rxn']:
            metric_optima.append(entry['gamma_new'])
            data_optimum.append(entry)

    label_colors = {'all': 'orange', 'observed':'green', 'optimum':'lightskyblue', 'ICSD':'forestgreen'}

    for data, label in zip([data_true, data_optimum], ['observed', 'optimum']):
    # for data, label in zip([data_true, data_optimum], ['True', 'Optimum']):
        sorted_data = np.sort(data)
        color = label_colors[label]
        cdf = np.arange(1, len(sorted_data) + 1) / len(sorted_data)
        ax.plot(sorted_data, cdf, label=f'{label}', linewidth=4, color = color)

    ax.set_xlim(-0.3,0.3)
    ax.hlines(y=0.96, xmin = -0.3, xmax = 0.3, color = 'deepskyblue', linestyle = '--', linewidth = 2)
    ax.hlines(y=0.92, xmin = -0.3, xmax = 0.3, color = 'forestgreen', linestyle = '--', linewidth = 2)

    ax.vlines(x=0.1, ymin = -0.2, ymax = 1.1, color = 'blue', linestyle = '--', linewidth = 3)
    ax.set_ylim(-0.01,1.01)
    ax.set_xlabel('$\Gamma$ (eV/atom)')
    ax.set_ylabel('Cumulative probability')
    ax.legend(loc = 'best', frameon = True, framealpha = 1)
    # ax.grid()
    ax.set_yticks([0,0.5,1.0])
    ax.set_yticklabels(['-0.0','0.5','1.0'])

    ax.set_xticks([-0.2,0,0.2])
    ax.tick_params(axis = 'y', right = False, length = 10, width = lw)
    ax.tick_params(axis = 'x', top = False, length = 10, width = lw)


def metric_hexbin_parity(ax, data_tm):
    """
    helper function for plotting Figure 2
    """

    # sns.set_theme(style="ticks", rc=set_rc_params())
    m_opt = []
    m_true = []
    for key in data_tm:
        for temperature in data_tm[key]:
            temp_data = data_tm[key][temperature]
            if temp_data:
                all_temp = [entry['gamma_new'] for entry in temp_data]
                true = [entry['gamma_new'] for entry in temp_data if entry['true_rxn']]
                opt = min(all_temp)
                m_opt.extend([opt]*len(true))
                m_true.extend(true)
    m_opt = np.asarray(m_opt).ravel()
    m_true = np.asarray(m_true).ravel()

    mask = m_true >= m_opt
    m_opt = m_opt[mask]
    m_true = m_true[mask]

    gridsize = (100,80)
    vmax = 60
    # Create Hexbin Parity Plot
    hb = ax.hexbin(m_opt, m_true, gridsize=gridsize, cmap='summer_r', mincnt=1,vmin=0, vmax = vmax)

    # Add Colorbar

    ax.plot([-10, 10], [-10, 10], linestyle = '--', color = 'black', lw=4, label="Parity Line")

    ax.set_xlim(-0.3,0.3)
    ax.set_ylim(-0.3,0.3)
    clip_path = Path([
        (-0.3, -0.3),
        (0.3, 0.3),
        (0.3, 0.3),
        (-0.3, 0.3),
        (-0.3, -0.3)
    ])
    hb.set_clip_path(PathPatch(clip_path, transform=ax.transData))

    
    # Labels and Title
    ax.set_xlabel('$\Gamma_{opt}$ (eV/atom)')
    ax.set_ylabel('$\Gamma_{obs}$ (eV/atom)')
    ax.set_yticks([-0.2,0,0.2])
    ax.set_xticks([-0.2,0,0.2])
    ax.tick_params(axis = 'y', right = False, length = 10, width = lw)
    ax.tick_params(axis = 'x', top = False, length = 10, width = lw)
    # ax.grid()
    return hb

# Figure 4 helpers

def binned_fraction_overlay(ax, x, scores, xedges,
                            color='cornflowerblue',
                            xlim=None, extend=True):
    """
    Markers at bin centers, line extended to plot limits.
    """
    x = np.asarray(x)
    scores = np.asarray(scores)
    labels = scores > 0.5

    # --- Fractions per bin ---
    frac_true = []
    for left, right in zip(xedges[:-1], xedges[1:]):
        mask = (x >= left) & (x < right)
        if mask.sum() > 0:
            frac_true.append(labels[mask].mean())
        else:
            frac_true.append(np.nan)

    frac_true = np.asarray(frac_true)

    # --- Bin centers ---
    centers = 0.5 * (xedges[:-1] + xedges[1:])

    # --- Line x/y (extended) ---
    x_line = centers
    y_line = frac_true

    if extend and xlim is not None:
        x_line = np.concatenate(([xlim[0]], centers, [xlim[1]]))
        y_line = np.concatenate(([frac_true[0]], frac_true, [frac_true[-1]]))

    # --- Plot line ---
    ax.plot(
        x_line, y_line,
        color=color, linewidth=5, zorder=3
    )

    # --- Plot markers only at centers ---
    ax.scatter(
        centers, frac_true,
        edgecolor=color, facecolor='white',
        s=100, linewidth=1.5, zorder=4
    )

    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticklabels(['', '', '50%', '', '100%'])

    if xlim is not None:
        ax.set_xlim(xlim)

def make_heatmap(ax, x, y, title, xlabel, ylabel,
                 xlim, ylim, bins, vmax, 
                 return_bins=False):
    """
    Create a heatmap with white squares for zero counts.
    """

    ax.set_box_aspect(1)

    heatmap, xedges, yedges = np.histogram2d(
        x, y, bins=bins, range=[xlim, ylim]
    )

    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    heatmap_masked = np.ma.masked_where(heatmap == 0, heatmap)

    cmap = plt.cm.summer_r.copy()
    cmap.set_bad(color='white')

    im = ax.imshow(
        heatmap_masked.T,
        cmap=cmap,
        origin="lower",
        extent=extent,
        aspect="auto",
        vmin=1,
        vmax=vmax
    )

    ax.set_title(title)
    ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)

    ax.set_yticks([0, 0.5, 1])
    ax.set_xticks([0, 0.1, 0.2, 0.3, 0.4])
    ax.set_yticklabels(['', '0.5', '1.0'])
    ax.set_xticklabels(['0.0', '0.1', '0.2', '0.3', '0.4'])

    if return_bins:
        return im, xedges
    return im


def get_data(gamma_gen, key = 'hull_energy', lim = None, models = ['PU-CGNF', 'SynthNN','PU-CGCNN','SynCoTrain']):

    keys = [k for k in gamma_gen]
    if key == 'gamma':
        if lim: 
            keys = [k for k in gamma_gen if gamma_gen[k]['hull_energy']<lim]
    y = [[gamma_gen[k]['predictions'][model] for k in keys] for model in models]
    x = [gamma_gen[k][key] for k in keys]
    return x, y


# Figure 6 helpers
def model_hist(ax, data, model_name,threshold=0.5, bin_number = 25, sz = sz):
    # Define consistent bins — fixed width and full coverage of data
    bins = np.linspace(0, 1, bin_number)  # 30 uniform-width bins

    # Plot the histogram, returning patches
    counts, edges, patches = ax.hist(data, bins=bins, edgecolor='black', alpha=0.8, linewidth = 0.8*lw)

    # Loop over bins and recolor based on threshold
    for patch, left, right in zip(patches, edges[:-1], edges[1:]):
        bin_center = (left + right) / 2
        if bin_center > threshold:
            patch.set_facecolor('steelblue')   # above threshold
        else:
            patch.set_facecolor('tomato')  # below threshold

    # Add labels or styling if you like
    ax.set_xlabel('Score')
    ax.set_xticks([0,0.5,1])
    ax.set_xticklabels(['0','0.5','1'])
    ax.set_ylabel('No. of materials')
    ax.set_xlim(0,1)
    if model_name == 'SynthNN':
        ax.set_ylim(0,1500)
        ax.text(0.06,1425, model_name, size = 0.9*sz)
    if model_name == 'PU-CGNF':
        ax.set_ylim(0,480)
        ax.text(0.025,430, model_name, size = 0.9*sz)
    if model_name == 'PU-CGCNN':
        ax.set_ylim(0,480)
        ax.text(0.025,430, model_name, size = 0.9*sz)
    if model_name == 'SynCoTrain':
        ax.set_ylim(0,480)
        ax.text(0.025,430, model_name, size = 0.9*sz)
    


def plot_stacked_hist(gamma_gen, ax, cmap_name='viridis'):
    g0 = [] 
    g1 = [] 
    g2 = [] 
    g3 = [] 
    g4 = [] 
    for key in gamma_gen: 
        pred_frac = 0 
        for model in ['SynthNN','PU-CGNF','PU-CGCNN','SynCoTrain']: 
            if abs(gamma_gen[key]['predictions'][model]>0.5): 
                pred_frac += 1 
        if pred_frac == 0: 
            g0.append(gamma_gen[key]['hull_energy']) 
        if pred_frac == 1: 
            g1.append(gamma_gen[key]['hull_energy']) 
        if pred_frac == 2: 
            g2.append(gamma_gen[key]['hull_energy']) 
        if pred_frac == 3:     
            g3.append(gamma_gen[key]['hull_energy']) 
        if pred_frac == 4: 
            g4.append(gamma_gen[key]['hull_energy'])
    cols = mpl.colormaps[cmap_name].resampled(20)
    colors = [cols(0.1), cols(0.3), cols(0.5), cols(0.7), cols(0.95)]

    # choose bins (you can adjust)
    all_data = np.concatenate([g0, g1, g2, g3, g4])
    bins = np.linspace(0, 0.5, 60)

    ax.hist(
        [g0, g1, g2, g3, g4],
        bins=bins,
        stacked=True,
        color=colors,
        label=['0/4', '1/4', '2/4', '3/4', '4/4'],
        edgecolor='black',
        linewidth=0.8*lw
    )
    ax.vlines(x = [np.mean(g0), np.mean(g1), np.mean(g2), np.mean(g3), np.mean(g4)], colors = colors, linestyle = '--', ymin = 0, ymax = 170, lw = lw)
    ax.set_ylim(0,170)

    ax.set_xlabel("$\\it{E}$$_{hull}$ (eV/atom)")
    ax.set_ylabel("No. of materials")
    ax.set_xlim(0, 0.5)
    handles, labels = ax.get_legend_handles_labels()

# Reorder however you want — e.g., reverse:
    handles = handles[::-1]
    labels = labels[::-1]

    ax.legend(handles, labels, frameon = True)