import os
import pandas as pd
import numpy as np
from scipy.stats import spearmanr, pearsonr
from scipy.stats import normaltest, ttest_ind
from statannotations.Annotator import Annotator

#%%  Plotting Defaults

import textwrap
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.ticker
import matplotlib.colors
import seaborn as sns
sns.set_theme(palette='tab10', style='white')

mpl.rcParams['figure.dpi'] = 100
mpl.rcParams['font.size'] = 18
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams["errorbar.capsize"] = 3.5
#plt.style.use('ggplot')

def normalize_by_batch(dfo):

    #dfo = dfo[dfo['batch'] != 'B4']

    index_keys = dfo.index.values
    if isinstance(index_keys[0], tuple):
        df_keys = set([x[0] for x in dfo.index.values])
        df_list = []
        for k in df_keys:
            df_list.append(dfo.loc[k])
    else:
        df_list = [dfo]

    df_list_norm = []
    for df in df_list:
        df_norm = pd.DataFrame(columns=dfo.columns)
        for b in df['batch'].unique():
            df_batch = df.loc[df['batch'] == b]
            if b == 'B4':
                b = 'B3'
            df_control = df.loc[(df['group'] == 'control') & (df['batch'] == b)].iloc[:, 4:]
            df_control = df_control.describe().iloc[1,:]

            df_batch_norm = df_batch.iloc[:, 4:] / df_control.values
            df_batch_norm.insert(0, 'subject', df_batch['subject'])
            df_batch_norm.insert(1, 'batch', df_batch['batch'])
            df_batch_norm.insert(2, 'group', df_batch['group'])
            df_batch_norm.insert(3, 'group2', df_batch['group2'])

            df_norm = df_norm.append(df_batch_norm)

        #df_norm.sort_index(inplace=True)
        df_norm.sort_values(by='group', inplace=True)
        df_list_norm.append(df_norm)
        del df_norm, df_batch

    if len(df_list) > 1:
        df_norm_all = pd.concat(df_list_norm, keys=df_keys)
    else:
        df_norm_all = df_list_norm[0]

    return df_norm_all

def boxplots(col, titlestr, df, savedir, batchnorm = False):
    # col       = string, column of the dataframe, df, to plot
    # titlestr  = string, tag to add to end of column name for plot title
    # df        = dataframe containing data to plot
    # savedir   = string, subdirectory under current path to save figures
    # batchnorm = Boolean, if True, apply batch normalization to flow data

    if batchnorm == True:
        df = normalize_by_batch(df)

    # Set-up subplot with dpi = 100
    px = 1 / plt.rcParams['figure.dpi']
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(1200 * px, 400 * px), dpi=100)
    #fig = plt.figure(figsize=(500*px, 500*px), dpi=100)

    df.sort_values(by='batch', inplace=True)
    titlestr = col + titlestr

    def create_boxplot(ax, my_order, xlabs, cat_col, col, df, cpalette, xft):
        sns.boxplot(ax=ax, x=cat_col, y=col, data=df, order=my_order, palette=cpalette, fliersize=10)
        ax.set_xticklabels(xlabs, fontsize=xft)
        ax.set_title(titlestr, fontsize=16)
        ax.set_xlabel('')

    # Plot all groups
    my_order = ['control','drug','fusD','fusS','fusD_drug','fusS_drug']
    xlabs = ['Control', 'αPD-1', 'fusD','fusS','fusD + αPD-1', 'fusS + αPD-1']
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = 'Paired'
    create_boxplot(ax0, my_order, xlabs, 'group', col, df, cpalette, 11)

    # Plot FUS vs Combination
    my_order = ['control', 'drug', 'fus', 'combination']
    xlabs = ['Control','αPD-1','fus', 'fus + αPD-1'] # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=6) for l in xlabs]
    cpalette = ['#a6cee3','#1f78b4','#33a02c','#e31a1c']
    create_boxplot(ax1, my_order, xlabs, 'group2', col, df, cpalette, 13)

    # Plot Sparse vs Dense
    my_order = ['control', 'drug', 'fusD', 'fusS']
    xlabs = ['Control','αPD-1','fusD', 'fusS'] # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3','#1f78b4','#33a02c','#e31a1c']
    create_boxplot(ax2, my_order, xlabs, 'group3', col, df, cpalette, 13)

    sp = os.path.join(os.getcwd(), savedir)
    try:
        os.mkdir(sp)
    except Exception:
        pass

    titlestr = titlestr.replace('%', 'Per')
    titlestr = titlestr.replace('/', ' to ')
    plt.savefig(sp + titlestr + '_box.png', dpi='figure')
    plt.show()

#%% Dot Plot Function

from pandas.api.types import CategoricalDtype
def dotplots(col, titlestr, df, savedir, batchnorm = False, dothue = 'tab:blue'):
    # col       = string, column of the dataframe, df, to plot
    # titlestr  = string, tag to add to end of column name for plot title
    # df        = dataframe containing data to plot
    # savedir   = string, subdirectory under current path to save figures
    # batchnorm = Boolean, if True, apply batch normalization to flow data
    # dothue    = string, single color, or 'by_cat', or 'batch', or some other df col
    # cpalette  = string, if dothue is not a single color, this is colormap used


    if batchnorm == True:
        df = normalize_by_batch(df)

    df.sort_values(by='batch', inplace=True)
    titlestr = col + titlestr

    # Set-up subplot with dpi = 100
    px = 1 / plt.rcParams['figure.dpi']
    fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(1200*px, 400*px), dpi=100)

    def create_dotplot(ax, my_order, xlabs, cat_col, col, df, dothue, cpalette, xft):

        cat_group_order = CategoricalDtype(my_order, ordered=True)
        df[cat_col] = df[cat_col].astype(cat_group_order)
        df.sort_values(cat_col, inplace=True)

        sns.pointplot(cat_col, col, data=df, ax=ax, dodge=True, join=False, ci='sd', errwidth=.75,
                      capsize=.25, markers='s', color='black', order=my_order)

        if dothue == 'tab:blue':
            sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, color='tab:blue', order=my_order)
        else:
            if (dothue == 'by_cat'):
                sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, hue=cat_col, order=my_order,
                              palette=cpalette)
                ax.get_legend().remove()
            elif (dothue == 'batch'):
                sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, hue='batch', order=my_order,
                              palette=cpalette)
                sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, -.3), ncol=3, title=None, frameon=False)


        ax.set_xticklabels(xlabs, fontsize=xft)
        ax.set_title(titlestr, fontsize=16)
        ax.set_xlabel('')

    # Plot all groups
    my_order = ['control', 'drug', 'fusD', 'fusS', 'fusD_drug', 'fusS_drug']
    xlabs = ['Control', 'αPD-1', 'fusD', 'fusS', 'fusD + αPD-1', 'fusS + αPD-1']
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = 'Paired'
    create_dotplot(ax0, my_order, xlabs, 'group', col, df, dothue, cpalette, 11)
    ax0.set_xlabel('')

    # Plot FUS vs Combination
    my_order = ['control', 'drug', 'fus', 'combination']
    xlabs = ['Control', 'αPD-1', 'fus', 'fus + αPD-1']  # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3', '#1f78b4', '#33a02c', '#e31a1c']
    create_dotplot(ax1, my_order, xlabs, 'group2', col, df, dothue, cpalette, 13)
    ax1.set_xlabel('')

    # Plot Dense vs Sparse
    my_order = ['control', 'drug', 'fusD', 'fusS']
    xlabs = ['Control', 'αPD-1', 'fusD', 'fusS']  # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3', '#1f78b4', '#33a02c', '#e31a1c']
    create_dotplot(ax2, my_order, xlabs, 'group3', col, df, dothue, cpalette, 13)
    ax2.set_xlabel('')

    sp = os.path.join(os.getcwd(), savedir)
    try:
        os.mkdir(sp)
    except Exception:
        pass

    titlestr = titlestr.replace('%', 'Per')
    titlestr = titlestr.replace('/', ' to ')
    plt.savefig(sp + titlestr + '_box.png', dpi='figure')
    plt.show()


#%% Combined Box and Dot Plots

def get_pvals(cat, col, df):
    grps = list(df.groupby(cat).indices.keys())
    dfp = df.pivot(index='subject', columns=cat, values=col)
    pairs = []
    pvals = []
    for j in range(len(grps)):
        for k in range(j + 1, len(grps)):
            st, p = ttest_ind(dfp[grps[j]], dfp[grps[k]], nan_policy='omit', alternative='two-sided')
            if p <= 0.05:
                pairs.append([grps[j], grps[k]])
                pvals.append(p)
    return pairs, pvals

def create_boxplot(ax, my_order, xlabs, cat_col, col, df, cpalette, xft, titlestr):
    plot_params = {
        'data': df,
        'x': cat_col,
        'y': col,
        'order': my_order,
        'palette': cpalette,
        'fliersize': 10
    }

    pairs, pvals = get_pvals(cat_col, col, df)

    sns.boxplot(ax=ax, **plot_params)
    ax.set_xticklabels(xlabs, fontsize=xft)
    ax.set_title(titlestr, fontsize=14.5)
    ax.set_xlabel('')
    if len(pairs) > 0:
        annotator = Annotator(ax, pairs, **plot_params)
        annotator.configure(line_offset=-2, text_offset=-7)
        annotator.set_pvalues(pvals)
        annotator.annotate()

def create_dotplot(ax, my_order, xlabs, cat_col, col, df, dothue, cpalette, xft, titlestr):
    plot_params = {
        'data': df,
        'x': cat_col,
        'y': col,
        'order': my_order,
    }

    pairs, pvals = get_pvals(cat_col, col, df)

    # sort dataframe by my_order of the cat_col variable
    cat_group_order = CategoricalDtype(my_order, ordered=True)
    df[cat_col] = df[cat_col].astype(cat_group_order)
    df.sort_values(cat_col, inplace=True)

    ydata = df.groupby(cat_col).mean()
    errdata = df.groupby(cat_col).std()
    ax.errorbar(x=my_order, y=ydata[col], yerr=errdata[col], fmt='s', color='black',
                ecolor='lightgray', elinewidth=4, capsize=0, markersize=10)

    #sns.pointplot(ax=ax, dodge=True, join=False, ci='sd', errwidth=.75,
    #              capsize=.25, markers='s', color='black', **plot_params)

    if dothue == 'tab:blue':
        sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=8, color='tab:blue')
    else:
        if (dothue == 'by_cat'):
            sns.stripplot(jitter=True, ax=ax, size=8, hue=cat_col, palette=cpalette, **plot_params)
            ax.get_legend().remove()
        elif (dothue == 'batch'):
            sns.stripplot(jitter=True, ax=ax, size=8, hue='batch', palette=cpalette, **plot_params)
            sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, -.3), ncol=3, title=None, frameon=False)


    ax.set_xticklabels(xlabs, fontsize=xft)
    ax.set_title(titlestr, fontsize=14.5)
    ax.set_xlabel('')
    if len(pairs) > 0:
        annotator = Annotator(ax, pairs, **plot_params)
        annotator.configure(line_offset=0, text_offset=-7)
        annotator.set_pvalues(pvals)
        annotator.annotate()

def boxdotplots(col, titlestr, df, savedir, batchnorm = False, dothue = 'tab:blue'):
    if batchnorm == True:
        df = normalize_by_batch(df)

    # Set-up subplot with dpi = 100
    px = 1 / plt.rcParams['figure.dpi']
    gs= {'wspace': .3, 'hspace': .25}
    fig, (axT, axB) = plt.subplots(nrows=2, ncols=3, figsize=(1200 * px, 900 * px), dpi=100, gridspec_kw=gs)
    #fig = plt.figure(figsize=(500*px, 500*px), dpi=100)

    df.sort_values(by='batch', inplace=True)
    titlestr = col + titlestr

    ## BOX PLOTS
    # Plot all groups
    my_order = ['control','drug','fusD','fusS','fusD_drug','fusS_drug']
    xlabs = ['Control', 'αPD-1', 'fusD','fusS','fusD + αPD-1', 'fusS + αPD-1']
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = 'Paired'
    create_boxplot(axT[0], my_order, xlabs, 'group', col, df, cpalette, 11, titlestr)
    ylims = axT[0].get_ylim()

    # Plot FUS vs Combination
    my_order = ['control', 'drug', 'fus', 'combination']
    xlabs = ['Control','αPD-1','fus', 'fus + αPD-1'] # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3','#1f78b4','#33a02c','#e31a1c']
    create_boxplot(axT[1], my_order, xlabs, 'group2', col, df, cpalette, 13, titlestr)

    # Plot Sparse vs Dense
    my_order = ['control', 'drug', 'fusD', 'fusS']
    xlabs = ['Control','αPD-1','fusD', 'fusS'] # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3','#1f78b4','#33a02c','#e31a1c']
    create_boxplot(axT[2], my_order, xlabs, 'group3', col, df, cpalette, 13, titlestr)

    ## DOT PLOTS
    # Plot all groups
    my_order = ['control', 'drug', 'fusD', 'fusS', 'fusD_drug', 'fusS_drug']
    xlabs = ['Control', 'αPD-1', 'fusD', 'fusS', 'fusD + αPD-1', 'fusS + αPD-1']
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = 'Paired'
    create_dotplot(axB[0], my_order, xlabs, 'group', col, df, dothue, cpalette, 11, titlestr)

    # Plot FUS vs Combination
    my_order = ['control', 'drug', 'fus', 'combination']
    xlabs = ['Control', 'αPD-1', 'fus', 'fus + αPD-1']  # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3', '#1f78b4', '#33a02c', '#e31a1c']
    create_dotplot(axB[1], my_order, xlabs, 'group2', col, df, dothue, cpalette, 13, titlestr)

    # Plot Dense vs Sparse
    my_order = ['control', 'drug', 'fusD', 'fusS']
    xlabs = ['Control', 'αPD-1', 'fusD', 'fusS']  # df.groupby('group').mean().index
    xlabs = [textwrap.fill(l, width=7) for l in xlabs]
    cpalette = ['#a6cee3', '#1f78b4', '#33a02c', '#e31a1c']
    create_dotplot(axB[2], my_order, xlabs, 'group3', col, df, dothue, cpalette, 13, titlestr)

    for ax in axT:
        ax.set(ylim=ylims)
    for ax in axB:
        ax.set(ylim=ylims)

    sp = os.path.join(os.getcwd(), savedir)
    try:
        os.mkdir(sp)
    except Exception:
        pass

    titlestr = titlestr.replace('%', 'Per')
    titlestr = titlestr.replace('/', ' to ')
    plt.savefig(sp + titlestr + '.png', dpi='figure')
    plt.show()

#%%  PairPlots with HeatMap
def pairplot_heat(mrtiColumns, flowColumns, datadf, titlestr, savedir):
    def annotateplt(x=None, y=None, data=None, color=None, label=None):
        ax = plt.gca()
        ymax = data[ax.get_ylabel()].max()
        ymin = data[ax.get_ylabel()].min()
        ax.set_ylim(ymin - ymax * .25, ymax * 1.25)

        data2 = data[[ax.get_xlabel(), ax.get_ylabel()]].dropna()
        Srho, Sp = spearmanr(data2[ax.get_xlabel()], data2[ax.get_ylabel()])
        Prho, Pp = pearsonr(data2[ax.get_xlabel()], data2[ax.get_ylabel()])

        if Sp <= 0.05:
            text = f"Sp = {Srho:.2f}*"
            fw = 'bold'
        else:
            text = f"Sp = {Srho:.2f}"
            fw = 'normal'
        ax.text(.3, 1.15, text, transform=ax.transAxes, fontsize=10, horizontalalignment='left',
                verticalalignment='center', fontweight=fw)

        if Pp <= 0.05:
            text = f"R = {Prho:.2f}*"
            fw = 'bold'
        else:
            text = f"R = {Prho:.2f}"
            fw = 'normal'
        ax.text(.3, 1.05, text, transform=ax.transAxes, fontsize=10, horizontalalignment='left',
                verticalalignment='center', fontweight=fw)

        normal = plt.Normalize(-1, 1)
        bgcol = plt.cm.RdBu_r(normal(Prho))
        ax.set_facecolor(bgcol)

        return ax

    g = sns.pairplot(datadf, x_vars=mrtiColumns, y_vars=flowColumns, height=2, aspect=1, kind='reg', hue='regcol',
                     palette=['#000000', '#000000'], dropna=True)
    g.map(annotateplt, data=datadf)
    plt.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
    plt.tight_layout()

    sp = os.path.join(os.getcwd(), savedir)
    try:
        os.mkdir(sp)
    except Exception:
        pass

    titlestr = titlestr.replace('%', 'Per')
    titlestr = titlestr.replace('/', ' to ')
    plt.savefig(sp + titlestr + '_scatterHeat.png', dpi='figure')
    plt.show()

#%% PairPlot by Group

def pairplot(mrtiColumns, flowColumns, datadf, titlestr, savedir, my_hue=None, cpalette=['#33a02c', '#e31a1c']):

    def annotateplt(x=None, y=None, hue=None, data=None, color=None, gr_vals=None, order=None, label=None, **kwargs):
        ax = plt.gca()
        ymax = data[ax.get_ylabel()].max()
        ymin = data[ax.get_ylabel()].min()
        ax.set_ylim(ymin-ymax*.25, ymax*1.25)

       # gr_vals = pd.unique(data[hue].values)
        x_locs = [.05, .6]


        for i in range(len(gr_vals)):
            print(gr_vals[i] + ':' + color[i])
            datagr = data.loc[data[hue] == gr_vals[i]]
            datagr2 = datagr[[ax.get_xlabel(), ax.get_ylabel()]].dropna()

            Srho, Sp = spearmanr(datagr2[ax.get_xlabel()], datagr2[ax.get_ylabel()])
            Prho, Pp = pearsonr(datagr2[ax.get_xlabel()], datagr2[ax.get_ylabel()])

            if Sp <= 0.05:
                text = f"Sp = {Srho:.2f}*"
                fw = 'bold'
            else:
                text = f"Sp = {Srho:.2f}"
                fw = 'normal'
            ax.text(x_locs[i], 1.15, text, transform=ax.transAxes, fontsize=10, horizontalalignment='left',
                    verticalalignment='center', fontweight=fw, color=color[i])

            if Pp <= 0.05:
                text = f"R = {Prho:.2f}*"
                fw = 'bold'
                print(ax.get_xlabel() + ' & ' + ax.get_ylabel() + ': ' + text + ', p-val = ' + f"{Pp:.5f}*")
            else:
                text = f"R = {Prho:.2f}"
                fw = 'normal'
            ax.text(x_locs[i], 1.05, text, transform=ax.transAxes, fontsize=10, horizontalalignment='left',
                    verticalalignment='center', fontweight=fw, color=color[i])

        return ax

    #cpalette = ['#33a02c', '#e31a1c']

    #datadf = datadf.sort_values(by=my_hue)
    #gr_order = np.unique(datadf[my_hue].values)

    # sort dataframe by my_order of the cat_col variable
    cat_group_order = CategoricalDtype(['fus', 'combination','fusD','fusS'], ordered=True)
    datadf[my_hue] = datadf[my_hue].astype(cat_group_order)
    datadf.sort_values(my_hue, inplace=True)
    gi = pd.unique(datadf[my_hue].values).astype(object)
    print('pairplot: ' + gi)

    g = sns.pairplot(datadf, x_vars=mrtiColumns, y_vars=flowColumns, height=2, aspect=1,
                     hue=my_hue, kind='reg', palette=cpalette, dropna=True)

    g.map(annotateplt, data=datadf, color=cpalette, hue=my_hue, gr_vals=gi)
   # g._legend.set_bbox_to_anchor((1, 0.9), borderaxespad=0)
    plt.legend(bbox_to_anchor=(1.02, 0.55), loc='upper left', borderaxespad=0)
    plt.tight_layout()

    sp = os.path.join(os.getcwd(), savedir)
    try:
        os.mkdir(sp)
    except Exception:
        pass

    titlestr = titlestr.replace('%', 'Per')
    titlestr = titlestr.replace('/', ' to ')
    plt.savefig(sp + titlestr + '_scatter_' + my_hue + '.png', dpi='figure')
    plt.show()



# def create_dotplot(ax, my_order, xlabs, cat_col, col, df, dothue, cpalette, xft):
#     cat_group_order = CategoricalDtype(my_order, ordered=True)
#     df[cat_col] = df[cat_col].astype(cat_group_order)
#     df.sort_values(cat_col, inplace=True)
#
#     sns.pointplot(cat_col, col, data=df, ax=ax, dodge=True, join=False, ci='sd', errwidth=.75,
#                   capsize=.25, markers='s', color='black', order=my_order)
#
#     if dothue == 'tab:blue':
#         sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, color='tab:blue', order=my_order)
#     else:
#         if (dothue == 'by_cat'):
#             sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, hue=cat_col, order=my_order,
#                           palette=cpalette)
#             ax.get_legend().remove()
#         elif (dothue == 'batch'):
#             sns.stripplot(cat_col, col, data=df, jitter=True, ax=ax, size=7, hue='batch', order=my_order,
#                           palette=cpalette)
#             sns.move_legend(ax, "lower center", bbox_to_anchor=(.5, -.3), ncol=3, title=None, frameon=False)
#
#     ax.set_xticklabels(xlabs, fontsize=xft)
#     ax.set_title(titlestr, fontsize=16)
#     ax.set_xlabel('')