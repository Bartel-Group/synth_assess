import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import CSS4_COLORS
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pydmclab.utils.handy import read_json
from pydmclab.core.comp import CompTools
from pydmclab.plotting.utils import set_rc_params
from synth_assess.plotting.plot_helpers import plot_hull, plot_rolling_stats, metric_cdf, metric_hexbin_parity, binned_fraction_overlay, make_heatmap, get_data, model_hist, plot_stacked_hist
from synth_assess.data.load import tm_entries, mp_data, mp_data_with_theoretical, gen_data, tm_rxns_with_gamma

set_rc_params()
COLORS = CSS4_COLORS # or some other palette
models = ['PU-CGNF', 'SynthNN','PU-CGCNN','SynCoTrain']
lw = 3.5
sz = 35

def plot_fig_1(sz = sz, p_sz = 250, figsize = (20,18)):
    tm_list = tm_entries()
    entries = [entry for entry in tm_list if entry['E_d'] != None]
    x = []
    for entry in entries:
        if entry['true_rxn'] == True:
            if entry['E_d']<0:
                x.append(0)
            else:
                x.append(entry['E_d'])
    y =  [entry['gamma_new'] for entry in entries if entry['true_rxn'] == True]



    # dummy data
    plt.rcParams.update({
    'font.size': sz,
    'axes.titlesize': sz,
    'axes.labelsize': sz,
    'xtick.labelsize': sz,
    'ytick.labelsize': sz,
    'legend.fontsize': 0.8*sz,
    'xtick.major.width': lw,
    'ytick.major.width': lw,
    'axes.linewidth': lw
    })

    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(2, 2, height_ratios = [1.2,1], hspace = 0.32, wspace = 0.15)  # 2 rows, 2 columns, equal height

    # top subplot spanning both columns
    ax1 = fig.add_subplot(gs[0, :])
    plot_rolling_stats(ax1, x, y)
    ax1.scatter([0.0311],[0.102], color = 'red', s = p_sz)
    ax1.scatter([0.0325],[0.0407], color = 'royalblue', s = p_sz)
    ax1.text(x=-0.02, y= 0.3, s = 'a', fontweight = 'bold', fontsize = 45)
    ax1.tick_params(axis='y', right = False, length = 10, width = lw)
    ax1.tick_params(axis='x', length = 10, top = False, width = lw)
    # bottom left
    ax2 = fig.add_subplot(gs[1, 0])
    plot_hull(ax2, 'LaNb$_3$O$_9$', p_sz=p_sz)
    ax2.set_yticks([-0.2, -0.1,0])    
    ax2.text(x=-0.43, y= 0.06, s = 'b', fontweight = 'bold', fontsize = 45)
    ax2.tick_params(axis='y', right = False, length = 10, width = lw)
    ax2.tick_params(axis='x', length = 10, top = False, width = lw)
    precs = ['Nb$_2$O$_5$','La$_2$O$_3$']
    ax2.set_xlabel(xlabel= '$\it{x}$ in (1 - $\it{x})$'+precs[0] + ' + ' + '$\it{x}$' + precs[1],size = 0.9*sz)

    

    # ax2.set_title("Bottom left")

    # bottom right
    ax3 = fig.add_subplot(gs[1, 1])
    plot_hull(ax3, 'LaGaO$_3$', ylabel= '', p_sz = p_sz)
    ax3.set_yticks([-0.2, -0.1,0])
    ax3.set_yticklabels(['','',''])
    ax3.tick_params(axis='y', right = False, length = 10, width = lw)
    ax3.tick_params(axis='x', length = 10, top = False, width = lw)
    precs = ['Ga$_2$O$_3$','La$_2$O$_3$']
    ax3.set_xlabel(xlabel= '$\it{x}$ in (1 - $\it{x})$'+precs[0] + ' + ' + '$\it{x}$' + precs[1],size = 0.9*sz)


    # ax3.set_title("Bottom right")

    plt.tight_layout()
    plt.show()



def plot_fig_2(sz = sz):
    tm_list = tm_entries()
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.8*sz,
        'xtick.major.width': lw,
        'ytick.major.width': lw,
        'axes.linewidth': lw
    })  
    fig, axes = plt.subplots(1, 2, figsize=(16, 6.375), constrained_layout = True)
    hb = metric_hexbin_parity(axes[0], data_tm)
    axes[0].text(x = -0.58,y = 0.3, s = 'a', fontweight = 'bold', fontsize = 45)
    metric_cdf(axes[1], tm_list)
    axes[0].tick_params(axis='y', right = False, length = 10, width = lw)
    axes[0].tick_params(axis='x', length = 10, top = False, width = lw)
    axes[1].text(x = -0.65,y = 1.0, s = 'b', fontweight = 'bold', fontsize = 45)
    axes[1].tick_params(axis='y', right = False, length = 10, width = lw)
    axes[1].tick_params(axis='x', length = 10, top = False, width = lw)
    ticks = np.linspace(0, hb.get_array().max(), 3)
    print(ticks)
    fig.colorbar(mappable= hb, ax=axes[0], cmap = 'summer_r', ticks = [0, 20, 40, 60], label = 'Reactions')


    plt.show()
    return


def plot_fig_3(sz = sz, xlim = [-0.1, 0.5], ylim = [0,175], ylim2 = [0,50], binwidth = 0.01, lim_color = 'blue'):
    gamma_gen = gen_data()
    x = [gamma_gen[k]['hull_energy'] for k in gamma_gen]
    x = [entry for entry in x if xlim[0]<= entry <= xlim[1]]
    x2 = [gamma_gen[k]['gamma'] for k in gamma_gen if gamma_gen[k]['hull_energy']<0.1]
    x2 = [entry for entry in x2 if xlim[0]<= entry <= xlim[1]]
    plt.rcParams.update({
            'font.size': sz,
            'axes.titlesize': sz,
            'axes.labelsize': sz,
            'xtick.labelsize': sz,
            'ytick.labelsize': sz,
            'legend.fontsize': sz,
            'xtick.major.width': lw,
            'ytick.major.width': lw,
            'axes.linewidth': lw
        })    
    fig, axes = plt.subplots(1,2, figsize = (16, 8), tight_layout=True)
    # bins = int((xlim[1] - xlim[0])/binwidth)
    bin_edges = np.arange(xlim[0], xlim[1] + binwidth, binwidth)
    # split left-hand data at E_hull = 0.1
    x_low  = [v for v in x if v < 0.1]
    x_high = [v for v in x if v >= 0.1]

    # plot >= 0.1 first (background)
    axes[0].hist(
        x_high,
        bins=bin_edges,
        edgecolor='black',
        color='lightgreen',
        alpha=0.8,
        linewidth=0.55*lw
    )

    # plot < 0.1 on top, same color as right-hand plot
    axes[0].hist(
        x_low,
        bins=bin_edges,
        edgecolor='black',
        color='green',        # same as right panel
        alpha=0.8,
        linewidth=0.55*lw
    )

    axes[0].set_xlim(xlim[0],xlim[1])
    axes[0].set_ylim(ylim[0],ylim[1])
    axes[0].set_ylabel('No. of materials')
    

    axes[0].set_xlabel('$\it{E}$$_{hull}$ (eV/atom)')
    axes[0].text(0.11,150,'$\it{E}$$_{hull}$ = 0.1 ' + r'$\frac{eV}{atom}$' , color = lim_color, size = int(sz*0.9))
    # axes[0].set_yticklabels(['','','',''])
    axes[0].vlines(x=0.1,ymin = 0, ymax = 500, color = lim_color, linewidth = lw, linestyle = '--')
    axes[0].tick_params(axis = 'y', right = False, length = 10, width = lw)
    axes[0].tick_params(axis = 'x', top = False, length = 10, width = lw)
    axes[0].set_yticks([0, 50, 100, 150])
    axes[0].set_yticklabels(['0','50','100','150'])
    axes[0].set_xticks([-0.1, 0, 0.1, 0.2, 0.3,0.4, 0.5])
    axes[0].set_xticklabels(['','0','','0.2','','0.4',''])
    axes[1].hist(x2, bins = bin_edges, edgecolor = 'black', color = 'green', alpha = 0.8, linewidth = 0.55*lw)
    axes[1].set_xlim(xlim[0],xlim[1])
    axes[1].set_ylim(ylim2[0],ylim2[1])
    axes[1].text(0.11,50,'$\Gamma_{opt}$ = 0.1 ' + r'$\frac{eV}{atom}$', color = lim_color, size = int(sz*0.9))
    axes[0].text(-0.33,175, s = 'a', fontweight = 'bold', fontsize = 45)
    axes[1].set_ylabel('No. of materials')
    # axes[1].set_yticklabels(['0','20','40','060'])
    axes[1].tick_params(axis = 'y', right = False, length = 10, width = lw)
    axes[1].tick_params(axis = 'x', top = False, length = 10, width = lw)
    axes[1].set_yticks([0,20,40,60])
    axes[1].set_yticklabels(['0','20','40','060'])
    axes[1].set_xticks([-0.1, 0, 0.1, 0.2, 0.3,0.4, 0.5])
    axes[1].set_xticklabels(['','0','','0.2','','0.4',''])
    axes[1].text(-0.33,60, s = 'b', fontweight = 'bold', fontsize = 45)
    axes[1].set_xlabel('$\Gamma_{opt}$ (eV/atom)')
    axes[1].vlines(x=0.1,ymin = 0, ymax = 500, color = lim_color, linewidth = 0.8*lw, linestyle = '--')



def plot_fig_4(key='hull_energy',
               models=['PU-CGNF', 'SynthNN', 'PU-CGCNN', 'SynCoTrain'],
               bins=12, vmax=100, score_line=None,
               ylim=(0, 1.0), lim=None, sz=sz):
    """
    Create 2x2 grid of heatmaps with binned fraction overlays.
    - One overlay point per heatmap x-bin
    - Points centered on each heatmap box
    """
    gamma_gen = gen_data()
    if key == 'hull_energy':
        xlim = (-0.0, 0.35)
    else:
        xlim=(-0.05, 0.35)
    # --- Data ---
    x_data, y_data = get_data(gamma_gen, key, lim, models=models)

    # --- Style ---
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.8 * sz,
        'axes.linewidth': lw
    })

    # --- Figure layout ---
    fig = plt.figure(figsize=(20, 21))
    spec = gridspec.GridSpec(2, 2, wspace=0.15, hspace=0.3)
    plots = []

    for i, model in enumerate(models):
        row, col = divmod(i, 2)
        ax_heat = fig.add_subplot(spec[row, col])

        x = x_data
        y = y_data[i]

        xl = '$\Gamma_{opt}$ (eV/atom)' if key == 'gamma' else r'$\mathit{E}_{hull}$ (eV/atom)'
        vmax_use = 50 if key == 'gamma' else vmax

        # --- Heatmap (return x-bin edges) ---
        im, xedges = make_heatmap(
            ax_heat, x, y,
            title='',
            xlabel=xl if row == 1 else '',
            ylabel=None,
            xlim=xlim,
            ylim=ylim,
            bins=bins,
            vmax=vmax_use,
            sz=sz,
            score_line=score_line,
            return_bins=True
        )
        plots.append(im)

        # Remove default left axis
        ax_heat.tick_params(left=False, labelleft=False)
        ax_heat.spines['left'].set_visible(False)
        ax_heat.spines['right'].set_visible(True)

        # --- Binned fraction overlay ---
        ax_frac = ax_heat.twinx()

        binned_fraction_overlay(
            ax_frac, x, y, xedges,
            color='blue',
            xlim=xlim
        )

        ax_frac.set_ylim(0, 1)
        ax_frac.tick_params(axis='y', right=False, labelright=False)

        ax_frac.spines['left'].set_visible(True)
        ax_frac.spines['left'].set_color('blue')
        ax_frac.tick_params(axis='y', colors='blue', width=lw,
                            labelleft=True, length=10)
        ax_frac.yaxis.set_ticks_position('left')

        if col == 0:
            ax_frac.set_ylabel('% pred. synthesizable', color='blue')
            ax_frac.yaxis.set_label_position('left')
            ax_heat.tick_params(right=True, labelright=False)

            if key != 'gamma':
                ax_frac.text(-0.06, 0.0, '0%', color='blue')
            else:
                ax_frac.text(-0.12, -0.02, '0%', color='blue')

        else:
            ax_frac.set_yticklabels([''] * 5)
            ax_frac.set_ylabel('')
            ax_heat.tick_params(axis='y', labelright=True,
                                right=True, colors='black', width=3)
            ax_heat.spines['right'].set_color('black')
            ax_heat.yaxis.set_label_position("right")
            ax_heat.yaxis.set_ticks_position("right")
            ax_heat.set_ylabel("Score", color='black')
            ax_heat.text(0.36, 0.0, '0.0')

        # --- Titles ---
        if key != 'gamma':
            xpos = 0.21 if model in ['PU-CGNF', 'SynthNN'] else 0.17
        else:
            xpos = 0.19 if model in ['PU-CGNF', 'SynthNN'] else 0.15

        ax_heat.text(xpos, 0.90, model)

        if row != 1:
            ax_heat.set_xlabel('')

        ax_heat.tick_params(axis='both', length=10, width=lw)
        ax_frac.tick_params(axis='y', length=10, width=lw)
        ax_heat.tick_params(axis='x', top=False, length=10)
        ax_frac.tick_params(axis='x', top=False, length=10)

    # --- Horizontal colorbar centered above top row ---
    ax_tl = fig.axes[0]
    ax_tr = fig.axes[1]

    pos_l = ax_tl.get_position()
    pos_r = ax_tr.get_position()

    cb_width = 1.2 * (pos_r.x1 - pos_l.x0)
    cb_height = 0.02
    cb_left = pos_l.x0 + (pos_r.x1 - pos_l.x0 - cb_width) / 2 + 0.2
    cb_bottom = pos_l.y1 + 0.05

    cbar_ax = fig.add_axes([cb_left, cb_bottom, cb_width, cb_height])

    ticks = [0, 10, 20, 30, 40, 50] if key == 'gamma' else [0, 20, 40, 60, 80, 100]

    cbar = fig.colorbar(
        plots[0],
        cax=cbar_ax,
        orientation='horizontal',
        ticks=ticks
    )

    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')
    cbar.ax.set_xticklabels(ticks[:-1] + [f"$\\geq${vmax_use}"])
    cbar.set_label("Number of targets", labelpad=4)

    return fig


## Fig 3 fns



def plot_fig_6(models=['PU-CGNF', 'SynthNN', 'PU-CGCNN', 'SynCoTrain'], cmap_name='viridis', sz = sz):
    gamma_gen = gen_data()
    fig = plt.figure(figsize=(20, 24))
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.9 * sz,
        'axes.linewidth': lw
    })

    # Outer GridSpec: 2 vertical regions (top block of 2 rows, bottom block of 1 row)
    outer_gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1], hspace=0.26, figure=fig)

    # Inner GridSpec for the top block (first 2 rows)
    top_gs = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer_gs[0], hspace=0.22, wspace=0.25)

    # Inner GridSpec for the bottom block (third row spanning both columns)
    bottom_gs = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec=outer_gs[1])

    # Define subplots
    ax1 = fig.add_subplot(top_gs[0, 0])
    ax2 = fig.add_subplot(top_gs[0, 1])
    ax3 = fig.add_subplot(top_gs[1, 0])
    ax4 = fig.add_subplot(top_gs[1, 1])
    ax5 = fig.add_subplot(bottom_gs[0])

    # Collect all axes
    axes = [ax1, ax2, ax3, ax4, ax5]

    # Example data & plotting loop
    x, model_data = get_data(gamma_gen, 'hull_energy')
    for i, ax in enumerate(axes[:-1]):  # first four subplots
        data = model_data[i]
        model = models[i]
        # if model == 'SynthNN':
        #     sz = 0.8*sz
        model_hist(ax, data, model, sz = sz)
        if i < 2:
            ax.set_xlabel('')
        if i in [1, 3]:
            ax.set_ylabel('')
        ax.tick_params(width = lw, length = 10, right = False, top = False)
    ax1.text(-0.4, 480, 'a', fontweight = 'bold', size = 1.1*sz)
    ax1.set_yticks([0,200,400])
    ax1.set_yticklabels(['0','200','400'])
    ax1.set_ylim(0,480)
    ax2.set_yticks([0,600,1200])
    ax2.set_yticklabels(['0','600','1200'])
    ax2.set_ylim(0,1600)
    ax2.tick_params(axis = 'y', labelsize = 0.9*sz)
    ax3.set_yticks([0,200,400])
    ax3.set_yticklabels(['0','200','400'])
    ax3.set_ylim(0,480)
    ax4.set_yticks([0,200,400])
    ax4.set_yticklabels(['0','200','400'])
    ax4.set_ylim(0,480)

    ax5.text(-0.09, 160, 'b', fontweight = 'bold', size = 1.1*sz)
    # Bottom subplot (row 3)
    plot_stacked_hist(gamma_gen, ax5, cmap_name=cmap_name)
    ax5.set_yticklabels(['','50','100','150'])
    plt.rcParams.update({'legend.fontsize':0.8*sz})
    plt.legend(frameon = True)
    plt.tight_layout()

    plt.show()





def plot_fig_s1(gamma_gen, sz=45, lw=lw):
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.6*sz,
        'xtick.major.width': lw,
        'ytick.major.width': lw,
        'axes.linewidth': lw
    })  
    fig, ax = plt.subplots(figsize = (10,8))
    ehs = [gamma_gen[k]['hull_energy'] for k in gamma_gen]
    gs = [gamma_gen[k]['gamma'] for k in gamma_gen]
    plot_rolling_stats(ax, ehs, gs)
    ax.set_ylabel('$\Gamma_{opt}$ (eV/atom)')


def plot_fig_s2(gamma_gen, figsize=(20,18), lower=False, sz = 35):
    """
    Plot pairwise scatter plots for given models with axes shown only on:
      - leftmost plots (y-axis)
      - top row plots (x-axis, labeled at top)
    """
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.9 * sz,
        'axes.linewidth': lw
    })
    m_list = ['PU-CGNF', 'SynthNN', 'PU-CGCNN', 'SynCoTrain']
    a, models = get_data(gamma_gen)
    models_dict = {m_list[i]:models[i] for i in range(4)}
    model_names = list(models_dict.keys())
    n = len(model_names)
    
    fig, axes = plt.subplots(4, 4, figsize=figsize, sharex=False, sharey=False, gridspec_kw={'hspace':0.2,'wspace':0.2})

    for i in range(n):
        for j in range(n):
            ax = axes[i, j]

            # Diagonal blank
            if i == j:
                ax.axis("off")
                continue

            # Skip redundant half (upper or lower)
            if (not lower and j <= i) or (lower and j >= i):
                ax.axis("off")
                continue

            xname, yname = model_names[j], model_names[i]
            x = models_dict[xname]
            y = models_dict[yname]

            ax.scatter(x, y, s=20, alpha=0.6)

            # --- Axes visibility logic ---

            # Leftmost column: show y-axis + ylabel
            ax.set_xlim(0,1)
            ax.set_ylim(0,1)
            ax.set_xticks([0,0.5,1])
            ax.set_yticks([0,0.5,1])
            ax.hlines(xmin = 0, xmax = 1, y = 0.5, color = 'blue', linewidth = 3, zorder = 5, linestyle = '--')
            ax.vlines(ymin = 0, ymax = 1, x = 0.5, color = 'blue', linewidth = 3, zorder = 5, linestyle = '--')

            ax.tick_params(top = False, right = False, length = 10, width =3)
            if j == i+1:
                ax.set_ylabel(yname)
                ax.set_xlabel(xname)
                ax.set_xticklabels(['0','','1'])
                ax.set_yticklabels(['0','','1'])
            else:
            #     ax.set_yticks([])
                ax.set_ylabel("")
                ax.set_xticklabels(['','',''])
                ax.set_yticklabels(['','',''])
                # ax.spines["left"].set_visible(False)
            if i == 0:
                # ax.xaxis.tick_top()
                ax.xaxis.set_label_position("top")
                ax.set_xlabel(xname)
            else:
                ax.set_xlabel("")


    plt.tight_layout()
    return fig, axes




def plot_fig_s3(gamma_gen, ylabel="Score", bins=12, vmax=100, score_line=None,
                       xlim=(-0.05, 0.35), ylim=(0, 1.0), sz=40):
    """
    Create 2x2 grid of heatmaps with moving fraction overlays.
    - Left column: colored left y-axis (%), no right axis.
    - Right column: no left axis, visible right y-axis (heatmap scale, labeled on right).
    - Colorbar centered vertically between rows.
    """

    # --- Data ---
    ex_data, ey_data = get_data(gamma_gen, key = 'hull_energy', lim = None, models=['TSDNN'])
    gx_data, gy_data = get_data(gamma_gen, key='gamma', lim = 0.1, models=['TSDNN'])


    # --- Style ---
    plt.rcParams.update({
        'font.size': sz,
        'axes.titlesize': sz,
        'axes.labelsize': sz,
        'xtick.labelsize': sz,
        'ytick.labelsize': sz,
        'legend.fontsize': 0.8 * sz,
        'axes.linewidth': lw
    })

    # --- Figure layout ---
    fig = plt.figure(figsize=(15, 15.5))
    spec = gridspec.GridSpec(2, 2, wspace=0.06, hspace=0.18)
    plots = []

    # for i, model in enumerate(models):
    #     row, col = divmod(i, 2)

    ax_heat = fig.add_subplot(spec[0, 0])
    x = ex_data
    y = ey_data[0]
    xl = r'$\mathit{E}_{hull}$ (eV/atom)'
    vmax = 100
    # --- Base heatmap ---
    im, xedges = make_heatmap(ax_heat, x, y, '', xl, None,
                        xlim = (0,0.35), ylim = ylim, bins = bins, vmax = vmax, sz = sz, score_line = None, 
                        return_bins=True)
    plots.append(im)

        # Remove default left spine/labels
    ax_heat.tick_params(left=False, labelleft=False)
    ax_heat.spines['left'].set_visible(False)
    ax_heat.spines['right'].set_visible(True)


    # --- Moving fraction overlay ---
    ax_frac = ax_heat.twinx()
    # print(model)
    binned_fraction_overlay(ax_frac, x, y, xedges ,color='blue', xlim = (0,0.35))
    ax_frac.set_ylim(0, 1)
        
    # Hide both sides initially
    ax_frac.tick_params(axis='y', right=False, labelright=False)

    ax_frac.spines['left'].set_visible(True)
    ax_frac.spines['left'].set_color('blue')
    # ax_frac.spines['left'].set_linewidth(3)
    ax_frac.tick_params(axis='y', colors='blue', width=lw, labelleft=True, length = 10)
    ax_frac.yaxis.set_ticks_position('left')
    ax_frac.yaxis.set_label_position('left')
    ax_frac.set_ylabel('% pred. synthesizable', color='blue')
    ax_frac.text(x=-0.05, y=0, s= '0%', color ='blue')
            # Hide right axis for heatmap
    ax_heat.tick_params(right=True, labelright=False)
    ax_heat.spines['right'].set_visible(True)

    gax_heat = fig.add_subplot(spec[0, 1])
    gx = gx_data
    gy = gy_data[0]
    xl = r'$\Gamma_{opt}$ (eV/atom)'
    vmax = 100
    # --- Base heatmap ---
    im, gxedges = make_heatmap(gax_heat, gx, gy, '', xl, None,
                        xlim = (-0.05,0.35), ylim = ylim, bins = bins, vmax = vmax, sz = sz, score_line = None,
                        return_bins = True)
    plots.append(im)

        # Remove default left spine/labels
    gax_heat.tick_params(left=False, labelleft=False)
    gax_heat.spines['left'].set_visible(False)
    gax_heat.spines['right'].set_visible(True)


    # --- Moving fraction overlay ---
    gax_frac = gax_heat.twinx()
    # print(model)
    binned_fraction_overlay(gax_frac, gx, gy, gxedges, color='blue', xlim = (-0.05,0.35))
    gax_frac.set_ylim(0, 1)
        
    # Hide both sides initially
    gax_frac.tick_params(axis='y', right=False, labelright=False)

    gax_frac.spines['left'].set_visible(True)
    gax_frac.spines['left'].set_color('blue')
    # ax_frac.spines['left'].set_linewidth(3)
    gax_frac.tick_params(axis='y', colors='blue',labelleft=False)
    gax_frac.set_yticklabels(['','','','',''])
    gax_frac.yaxis.set_ticks_position('left')
    gax_frac.yaxis.set_label_position('left')
    gax_heat.yaxis.set_label_position('right')
    gax_frac.set_ylabel('')
            # Hide right axis for heatmap
    gax_heat.tick_params(right=True, labelright=True)
    gax_heat.spines['right'].set_visible(True)
    gax_heat.set_ylabel("Synthesizability score", color='black')
    gax_heat.set_yticklabels(['0.0','0.5','1.0'])
    


        # --- Horizontal colorbar at top spanning full width ---
    # Position the colorbar axis above the top row
    # [left, bottom, width, height] in figure coords
        # --- Horizontal colorbar centered above the top two plots ---
    # Get positions of the two top-row axes
    ax_tl = fig.axes[0]   # top-left
    ax_tr = fig.axes[1]   # top-right

    pos_l = ax_tl.get_position()
    pos_r = ax_tr.get_position()

    # Define a smaller width for the colorbar (e.g., 35% of top-row span)
    cb_width =  1.2*(pos_r.x1 - pos_l.x0)
    cb_height = 0.02

    # Center the colorbar horizontally between the top two plots
    cb_left = pos_l.x0 + (pos_r.x1 - pos_l.x0 - cb_width) / 2 + 0.2
    cb_bottom = pos_l.y1 + 0.05   # slightly above the top axes

    cbar_ax = fig.add_axes([cb_left, cb_bottom, cb_width, cb_height])

    ticks = [0, 20, 40, 60, 80, 100]



    cbar = fig.colorbar(
        plots[0],
        cax=cbar_ax,
        orientation='horizontal',
        ticks=ticks
    )

    # Put ticks + label on top
    cbar.ax.xaxis.set_ticks_position('top')
    cbar.ax.xaxis.set_label_position('top')

    cbar.ax.set_xticklabels(ticks[:-1] + [f"$\geq${vmax}"])
    cbar.set_label("Number of targets", labelpad=4)


    return fig


def main():
    data_mp = mp_data()
    data_mp_theoretical = mp_data_with_theoretical()
    data_mp ={CompTools(entry['formula']).clean: entry for entry in data_mp}
    data_tm = tm_rxns_with_gamma()
    entries = tm_entries()
    # data_gen = read_json(os.path.join(DATADIR_gen, 'gen_data_with_pred.json'))
    gamma_gen = gen_data()

    return data_mp, data_mp_theoretical, data_tm, stability, entries, gamma_gen

if __name__ == "__main__":
    data_mp, data_mp_theoretical, data_tm, stability, entries, gamma_gen = main()
