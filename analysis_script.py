import os
import pandas as pd
import numpy as np


def newgroups(df):
    cond = [(df.group.str.fullmatch('control')),
            (df.group.str.fullmatch('drug')),
            (df.group.str.fullmatch('fusS')),
            (df.group.str.fullmatch('fusD')),
            (df.group.str.contains('_drug'))]
    values = ['control','drug','fus','fus','combination']
    df['group2'] = np.select(cond, values)

    col = df.pop('group2')
    df.insert(3, col.name, col)

    cond = [(df.group.str.fullmatch('control')),
            (df.group.str.fullmatch('drug')),
            (df.group.str.contains('fusS')),
            (df.group.str.contains('fusD'))]
    values = ['control', 'drug', 'fusS', 'fusD']
    df['group3'] = np.select(cond, values)

    col = df.pop('group3')
    df.insert(4, col.name, col)
    return df

# import tumor volume CSV data
tumor_voldf = pd.read_pickle('tumor_volume_data.pkl')
tumor_voldf = tumor_voldf[['subject','Start Volume','End Volume','End % Change in Volume']]

# import h&e necrosis csv data
histdf = pd.read_pickle('h&e_necrosis_seg.pkl')

# import MRTI stats CSV data
datapath = '/Users/sarajohnson/Box Sync/Work Documents/XFUS Research/Immunotherapy Project/Experiments/'
mrtidf = pd.read_csv(datapath + 'MRTIstats.csv')
mrtidf.drop('tumortemp', axis=1, inplace=True)
mrtidf = mrtidf.set_axis(['subject','FUS group', 'MR volume','MR volumeBV','TD240','T60','Hyper','TD240_BV','T60_BV','Hyper_BV'], axis=1, inplace=False)
cond = [(mrtidf['FUS group'].str.fullmatch('Sp')),(mrtidf['FUS group'].str.fullmatch('D'))]
values = ['fusS', 'fusD']
mrtidf['FUS group'] = np.select(cond, values)
mrtidf.drop(mrtidf.loc[mrtidf['subject'] == 15].index, axis=0, inplace=True)


# import .csv as dataframe
datapath = '~/Box Sync/Flow Cytometry Data/Compiled Results/'
df = pd.read_csv(datapath + 'Tumor_counts.csv')
tumordata = df #tumor counts
df = pd.read_csv(datapath + 'Spleen_counts.csv')
spleendata = df #spleen counts

# in tumordata, set the column name to final pop
tumordata.rename(columns=lambda x: x[:-8].split('/')[-1], inplace=True)
tumordata.rename(columns={'Un': 'subject', '*': 'batch', '*T': 'group'}, inplace=True)
tumordata.drop(df.tail(2).index, inplace=True)
tumordata['subject'] = [int(x[1:-4]) for x in tumordata['subject']]
tumordata = tumordata.merge(tumor_voldf, on='subject')
tumordata = tumordata.merge(mrtidf, how ='outer', on='subject')
tumordata = tumordata.merge(histdf, how ='outer', on='subject')
tumordata.iloc[:, -6:] = tumordata.iloc[:, -6:].fillna(0)
tumordata.sort_values(by='group', inplace=True)
tumordata = newgroups(tumordata)
#tumordata = tumordata[tumordata.subject != 'M2.fcs']

# in spleendata, set the column name to final pop
spleendata.rename(columns=lambda x: x[:-8].split('/')[-1], inplace=True)
spleendata.rename(columns={'Un': 'subject', '*': 'batch', '*T': 'group'}, inplace=True)
spleendata.drop(df.tail(2).index, inplace=True)
spleendata.loc[(spleendata.subject == 'M26b.fcs'), 'subject'] = 'M26.fcs'
spleendata['subject'] = [int(x[1:-4]) for x in spleendata['subject']]
spleendata = spleendata.merge(tumor_voldf, on='subject')
spleendata = spleendata.merge(mrtidf, how ='outer', on='subject')
spleendata.iloc[:, -6:] = spleendata.iloc[:, -6:].fillna(0)
spleendata.sort_values(by='group', inplace=True)
spleendata = newgroups(spleendata)
#spleendata = spleendata[spleendata.subject != 'M2.fcs']

# combine spleen & tumor data. Create DFs for '% of Live' and '% of CD45'
ad = pd.concat([tumordata, spleendata], keys=['tumor', 'spleen'])
end1 = ad.columns.get_loc('Single Cells-1')
start1 = ad.columns.get_loc('Start Volume')
ad2 = pd.concat([ad.iloc[:, :end1], ad.iloc[:, start1:], ad.iloc[:, end1:start1-1]], axis=1)
ad = ad2

l1 = ad.columns.get_loc('Single Cells-1')
cLive = ad.columns.get_loc('Live')+1
cCD45 = ad.columns.get_loc('CD45+')+1
ad_perLive = pd.concat([ad.iloc[:,:l1], ad.iloc[:,cLive:].div(ad.Live, axis=0)], axis=1)
ad_perCD45 = pd.concat([ad.iloc[:,:l1], ad.iloc[:,cCD45:].div(ad['CD45+'], axis=0)], axis=1)
ad_other = ad.iloc[:,:l1]
    #['subject','batch','group','group2','Start Volume','End Volume','End % Change in Volume']]


del tumordata, spleendata, df, ad2

#%% Set-up  MRTI columns

# calculate additional ratios
ad_other['TD240:Hyper'] = ad['TD240']/ad['Hyper']
ad_other['T60:Hyper'] = ad['T60']/ad['Hyper']
ad_other['TD240:Hyper_BV'] = ad['TD240_BV']/ad['Hyper_BV']
ad_other['T60:Hyper_BV'] = ad['T60_BV']/ad['Hyper_BV']
mrtiR_col = ['TD240:Hyper','T60:Hyper', 'TD240:Hyper_BV', 'T60:Hyper_BV']
ad_other[mrtiR_col] = ad_other[mrtiR_col].fillna(0)

# create column lists for later plotting
mrti_col = list(ad_other.columns[ad_other.columns.get_loc('TD240'):ad_other.columns.get_loc('Hyper_BV')+1].values)
mrtiColumns = mrti_col + mrtiR_col
infoColumns = list(ad_other.columns[0:ad_other.columns.get_loc('Start Volume')].values)
volColumns = list(ad_other.columns[ad_other.columns.get_loc('Start Volume'):ad_other.columns.get_loc('MR volumeBV')+1].values)

mrti_sub2 = [x for x in mrtiColumns if '_BV' in x]
mrti_sub1 = [x for x in mrtiColumns if x not in mrti_sub2]
mrti_subs = [mrti_sub1, mrti_sub2]

# Set-up H&E columns
ad_other['HE_AllNecrosis'] = ad_other['HE_Necrosis1']+ad_other['HE_Necrosis2']
histColumns = ['HE_Necrosis1','HE_Negative','HE_Necrosis2','HE_Live','HE_AllNecrosis']

excludeColumns = mrtiColumns + infoColumns + volColumns + histColumns



#%%  Batch Normalization Function
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

#%%  Boxplot Function

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
from scipy.stats import normaltest, ttest_ind
from statannotations.Annotator import Annotator

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
from scipy.stats import spearmanr, pearsonr

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


#%% Create all Plots

def create_plots(boxplot_titlestr, pairplot_titlestr, datadf, savedir, norm, dotplothue, mrti_col, flow_col, bxplt=True, pairplt=True):

    if bxplt == True:
        for col in flow_col:
            boxdotplots(col, boxplot_titlestr, datadf, savedir, batchnorm=norm, dothue=dotplothue)

    if pairplt == True:
        df2 = df.loc[(df['group2'] == 'fus') | (df['group2'] == 'combination')]
        df2['regcol'] = 'k'
        pairplot_heat(mrti_col, flow_col, df2, pairplot_titlestr, savedir)
        pairplot(mrti_col, flow_col, df2, pairplot_titlestr, savedir, my_hue='group2')
        pairplot(mrti_col, flow_col, df2, pairplot_titlestr, savedir, my_hue='group3', cpalette=['#8E44AD', '#E67E22' ])

        del df2

    del datadf

#%% TUMOR GROWTH
ad_tg = ad_other.copy()
ad_tg['% Change in Volume'] = ad['End % Change in Volume']
ad_tg['Change in Volume (mm3)'] = ad['End Volume'] - ad['Start Volume']
ad_tg['End Volume'] = ad['End Volume']
ad_tg['Start Volume'] = ad['Start Volume']

savedir = 'Figures_New2/Tumor Growth/'
mrti_col = mrti_sub2
norm = False
dotplothue = 'by_cat'
flow_col = [x for x in ad_tg.columns.values if x not in excludeColumns]
flow_col = flow_col + ['End Volume', 'Start Volume']

# Tumor graphs
df = ad_tg.copy()
df = df.loc['tumor']
boxplot_titlestr = ''
pairplot_titlestr = 'Tumor Growth'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=False)

del ad_tg, df
#%% VIABILITY

ad_via = ad_other.copy()
ad_via['Viability'] = ad['Live']/ad['Single Cells-1']
ad_via['CD45+/EpCAM+ Ratio'] = ad['CD45+']/ad['EpCAM+']/100
ad_via['% EpCAM+ of Live'] = ad_perLive['EpCAM+']
ad_via['% CD45+ of Live'] = ad_perLive['CD45+']
# ad_via['Tumor Cell Viability'] = ad['EpCAM+']/ad['Single Cells-1']
# ad_via['Tumor Cell Viability 2'] = ad['EpCAM+']/ad['EpCAM+ Cells']
# ad_via['Leukocyte Viability'] = ad['CD45+']/ad['Single Cells-1']
# ad_via['Leukocyte Viability 2'] = ad['CD45+']/ad['CD45+ Cells']

savedir = 'Figures_New2/Viability/'
mrti_col = mrti_sub2
norm = False
dotplothue = 'by_cat'
flow_col = [x for x in ad_via.columns.values if x not in excludeColumns]

# Tumor graphs
df = ad_via.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
#df['Tumor Cell Viability'] = np.where(df['Tumor Cell Viability'] > .5, np.nan, df['Tumor Cell Viability'])
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'Viability in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

del ad_via, df

#%% HE Viability

ad_hist = ad_other.copy()
ad_hist['Flow Viability'] = ad['Live']/ad['Single Cells-1']

savedir = 'Figures_New2/HE_Viability/'
flow_col = histColumns
mrti_col = ['TD240_BV','T60_BV','Hyper_BV','Flow Viability' + 'End % Change in Volume']
norm = False
dotplothue = 'by_cat'

# Tumor graphs
df = ad_hist.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
df = df.loc[df['HE_Necrosis1'] <50]
#df['Tumor Cell Viability'] = np.where(df['Tumor Cell Viability'] > .5, np.nan, df['Tumor Cell Viability'])
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'H&E Viability in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=False)



#%% Lymphocytes
ad_til = ad_other.copy()
# CD3+ T cells
ad_til['% CD3+ of Live'] = ad_perLive['CD3+ Lymphocytes']
#ad_til['% CD3+ of CD45+'] = ad_perCD45['CD3+ Lymphocytes']
ad_til['% CD3+/EpCAM+ Ratio'] = ad['CD3+ Lymphocytes']/ad['EpCAM+']
# B Cells
ad_til['% B cells of Live'] = ad_perLive['B220+ CD11c- B Cells']
#ad_til['% B cells of CD45+'] = ad_perCD45['B220+ CD11c- B Cells']
# NK Cells
ad_til['% NK Cells of Live'] = ad_perLive['NKp46+ NK Cells']
#ad_til['% NK Cells of CD45+'] = ad_perCD45['NKp46+ NK Cells']
# total Lymphocytes
ad_til['% Lymphocytes of Live'] = ad_til['% CD3+ of Live'] + ad_til['% B cells of Live'] + ad_til['% NK Cells of Live']
#ad_til['% Lymphocytes of CD45+'] = ad_til['% CD3+ of CD45+'] + ad_til['% B cells of CD45+'] + ad_til['% NK Cells of CD45+']

mrti_col = mrti_sub2
norm = False
dotplothue = 'by_cat'

# ad_til = ad_other.copy()
# ad_til['% B cells of Live'] = ad_perLive['B220+ CD11c- B Cells']
# ad_til['% B cells of CD45+'] = ad_perCD45['B220+ CD11c- B Cells']
# ad_til = ad_til.loc[ad_til['batch'] != 'B3']

# TUMOR
savedir = 'Figures_New2/Lymphocytes_Tumor/'
flow_col = [x for x in ad_til.columns.values if x not in excludeColumns]
df = ad_til.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'Lymphocytes in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

#%%
# SPLEEEN
savedir = 'Figures_New2/Lymphocytes_Spleen/'
flow_col = [x for x in ad_til.columns.values if x not in excludeColumns]
flow_col = [x for x in flow_col if 'CD45+' not in x]
flow_col = [x for x in flow_col if 'EpCAM' not in x]
df = ad_til.copy()
df = df.loc['spleen']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Spleen'
pairplot_titlestr = 'Lymphocytes in Spleen'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

del ad_til, df

#%% T cells

# CD4+ T cells
ad_tc = ad_other.copy()
#ad_tc['% CD4+ of Live'] = ad_perLive['CD4+ T Cells']
#ad_tc['% CD4+ of CD45'] = ad_perCD45['CD4+ T Cells']
ad_tc['% CD4+ of CD3+'] = ad['CD4+ T Cells']/ad['CD3+ Lymphocytes']
# CD8+ T cells
#ad_tc['% CD8+ of Live'] = ad_perLive['CD8+ T Cells']
#ad_tc['% CD8+ of CD45'] = ad_perCD45['CD8+ T Cells']
ad_tc['% CD8+ of CD3+'] = ad['CD8+ T Cells']/ad['CD3+ Lymphocytes']
# CD8+/CD4+ Ratio
ad_tc['CD8+/CD4+'] = ad['CD8+ T Cells']/ad['CD4+ T Cells']/100
# Activation
#ad_tc['% CD62L+ of CD4+'] = ad['CD62L+ CD4+ T Cells']/ad['CD4+ T Cells']
#ad_tc['% CD62L+ of CD8+'] = ad['CD62L+ CD8+ T Cells']/ad['CD8+ T Cells']

mrti_col = [mrti_sub2[0], mrti_sub2[2]]
norm = False
dotplothue = 'by_cat'

# TUMOR
savedir = 'Figures_New2/Tcells_Tumor/'
flow_col = [x for x in ad_tc.columns.values if x not in excludeColumns]
df = ad_tc.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'T cells in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

#%%
#subset plot
df3 = df.loc[(df['group2'] == 'fus') | (df['group2'] == 'combination')]
df3['regcol'] = 'k'
flow_col = ['% CD4+ of CD3+', '% CD8+ of CD3+', 'CD8+/CD4+']
pairplot_heat(mrti_col, flow_col, df3, pairplot_titlestr + '_v2', savedir)

# SPLEEEN
savedir = 'Figures_New2/Tcells_Spleen/'
flow_col = [x for x in flow_col if 'CD45+' not in x]
flow_col = [x for x in flow_col if 'EpCAM' not in x]
df = ad_tc.copy()
df = df.loc['spleen']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Spleen'
pairplot_titlestr = 'T cells in Spleen'

flow_col = ['% CD4+ of CD3+', '% CD8+ of CD3+', 'CD8+/CD4+']
create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

df3 = df.loc[(df['group2'] == 'fus') | (df['group2'] == 'combination')]
df3['regcol'] = 'k'
flow_col = ['% CD4+ of CD3+', '% CD8+ of CD3+', 'CD8+/CD4+']
pairplot_heat(mrti_col, flow_col, df3, pairplot_titlestr + '_v2', savedir)

del ad_tc, df

#%% PD-1 Expression

ad_pd1 = ad_other.copy()
#ad_pd1['% PD-1+ of EpCAM+'] = ad['PD-1+ Tumor Cells']/ad['EpCAM+']
# CD8 T cells
ad_pd1['% PD-1+ of CD8+'] = ad['PD-1+ CD8+ T Cells']/ad['CD8+ T Cells']
#ad_pd1['% PD-1+ CD8+ of CD3+'] = ad['PD-1+ CD8+ T Cells']/ad['CD3+ Lymphocytes']
#ad_pd1['% PD-1+ CD8+ of Live'] = ad_perLive['PD-1+ CD8+ T Cells']
#ad_pd1['% PD-1+ CD8+ of CD45+'] = ad_perCD45['PD-1+ CD8+ T Cells']
# CD4 T cells
ad_pd1['% PD-1+ of CD4+'] = ad['PD-1+ CD4+ T Cells']/ad['CD4+ T Cells']
#ad_pd1['% PD-1+ CD4+ of CD3+'] = ad['PD-1+ CD4+ T Cells']/ad['CD3+ Lymphocytes']
#ad_pd1['% PD-1+ CD4+ of Live'] = ad_perLive['PD-1+ CD4+ T Cells']
#ad_pd1['% PD-1+ CD4+ of CD45+'] = ad_perCD45['PD-1+ CD4+ T Cells']
# B cells
ad_pd1['% PD-1+ of B Cells'] = ad['PD-1+ B Cells']/ad['B220+ CD11c- B Cells']
#ad_pd1['% PD-1+ B Cells of CD45+'] = ad_perCD45['PD-1+ B Cells']
#ad_pd1['% PD-1+ B Cells of Live'] = ad_perLive['PD-1+ B Cells']
# T cells
ad_pd1['% PD-1+ of T Cells'] = ad['PD-1+ CD3+ Lymphocytes']/ad['CD3+ Lymphocytes']
#ad_pd1['% PD-1+ T Cells of CD45+'] = ad_perCD45['PD-1+ CD3+ Lymphocytes']
#ad_pd1['% PD-1+ T Cells of Live'] = ad_perLive['PD-1+ CD3+ Lymphocytes']
# Leukocytes
ad_pd1['% PD-1+ of CD45+'] = ad_perCD45['PD-1+ CD45+']
#ad_pd1['% PD-1+ CD45+ of Live'] = ad_perLive['PD-1+ CD45+']
# TILs
ad_pd1['% PD-1+ of TILs'] = (ad['PD-1+ CD3+ Lymphocytes']+ad['PD-1+ B Cells'])/(ad['CD3+ Lymphocytes']+ad['B220+ CD11c- B Cells'])
#ad_pd1['% PD-1+ TILs of CD45+'] = (ad_perCD45['PD-1+ CD3+ Lymphocytes']+ad_perCD45['PD-1+ B Cells'])
#ad_pd1['% PD-1+ TILs of Live'] = (ad_perLive['PD-1+ CD3+ Lymphocytes']+ad_perLive['PD-1+ B Cells'])


mrti_col = mrti_sub2
norm = False
dotplothue = 'by_cat'


# TUMOR
savedir = 'Figures_New2/PD-1_Tumor/'
flow_col = [x for x in ad_pd1.columns.values if x not in excludeColumns]
df = ad_pd1.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'PD-1 in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=False)

#%%

# SPLEEN
savedir = 'Figures_New2/PD-1_Spleen/'
flow_col = [x for x in ad_pd1.columns.values if x not in excludeColumns]
flow_col = [x for x in flow_col if 'CD45+' not in x]
flow_col = [x for x in flow_col if 'EpCAM' not in x]
# Spleen graphs
df = ad_pd1.copy()
df = df.loc['spleen']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Spleen'
pairplot_titlestr = 'PD-1 in Spleen 2'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=False)

del ad_pd1, df

#%% DCs
ad_mye = ad_other.copy()
#ad_mye['DCs of Live'] = ad_perLive['CD11c+ DCs']
#ad_mye['DCs of CD45+'] = ad_perCD45['CD11c+ DCs']
# B220+ DCs
#ad_mye['% B220+ DCs of CD45+'] = ad_perCD45['B220+']
#ad_mye['% B220+ DCs of Live'] = ad_perLive['B220+']
#ad_mye['% pDCs of Live'] = ad_perLive['pDCs (MHCII_int)']
#B220- DCs
#ad_mye['% B220- DCs of CD45+'] = ad_perCD45['B220-']
#ad_mye['% B220- DCs of Live'] = ad_perLive['B220-']
ad_mye['% cDCs of Live'] = ad_perLive['MHCII+ cDCs']
ad_mye['% cDCs of CD45+'] = ad_perCD45['MHCII+ cDCs']
#ad_mye['CD1/CD2 ratio'] = ad['cDC1s (CD8a+)']/ad['cDC2s (CD4+)']/100
ad_mye['% MHCII_hi of cDCs'] = ad['MHCII_hi, CD11c_hi']/ad['MHCII+ cDCs']

mrti_col = mrti_sub2
mrti_col = ['TD240_BV','Hyper_BV','TD240:Hyper_BV']
norm = False
dotplothue = 'by_cat'


# TUMOR
savedir = 'Figures_New2/DCs_Tumor/'
flow_col = [x for x in ad_mye.columns.values if x not in excludeColumns]
df = ad_mye.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'DCs in Tumor2'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=False)

#%%
# SPLEEN
savedir = 'Figures_New2/DCs_Spleen/'
flow_col = [x for x in ad_mye.columns.values if x not in excludeColumns]
flow_col = [x for x in flow_col if 'CD45+' not in x]
flow_col = [x for x in flow_col if 'EpCAM' not in x]
# Spleen graphs
df = ad_mye.copy()
df = df.loc['spleen']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Spleen'
pairplot_titlestr = 'DCs in Spleen'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=True)

df3 = df.loc[(df['group2'] == 'fus') | (df['group2'] == 'combination')]
df3['regcol'] = 'k'
df3['Neutrophils of Live'] = np.where(df3['Neutrophils of Live'] > 15, np.nan, df3['Neutrophils of Live'])
flow_col = ['Macs of Live', 'DCs of Live','Neutrophils of Live']
pairplot_heat(mrti_col, flow_col, df3, pairplot_titlestr + '_v2', savedir)

del ad_mye, df

#%%
ad_mye = ad_other.copy()
#ad_mye['Macs of Live'] = ad_perLive['F4_80+ Macs']
ad_mye['Macs of CD45+'] = ad_perCD45['F4_80+ Macs']
#ad_mye['Macs of CD11b+ Monocytes'] = ad['F4_80+ Macs']/ad['CD11b+ monocytes (Ly6G-)']
#ad_mye['Monocytes of Live'] = ad_perLive['Ly6C_hi, Ly6G- Inflam Monocyte (MDSC)']
ad_mye['Monocytes of CD45+'] = ad_perCD45['Ly6C_hi, Ly6G- Inflam Monocyte (MDSC)']
#ad_mye['Neutrophils of Live'] = ad_perLive['Ly6G+ Neutrophils']
ad_mye['Neutrophils of CD45+'] = ad_perCD45['Ly6G+ Neutrophils']

mrti_col = mrti_sub2
mrti_col = ['TD240_BV','Hyper_BV','TD240:Hyper_BV']
norm = False
dotplothue = 'by_cat'


# TUMOR
savedir = 'Figures_New2/Myeloid_Tumor/'
flow_col = [x for x in ad_mye.columns.values if x not in excludeColumns]
df = ad_mye.copy()
df = df.loc['tumor']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Tumor'
pairplot_titlestr = 'Myeloid in Tumor'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=False, pairplt=True)

#%%
# SPLEEN
savedir = 'Figures_New2/Myeloid_Spleen/'
flow_col = [x for x in ad_mye.columns.values if x not in excludeColumns]
flow_col = [x for x in flow_col if 'CD45+' not in x]
flow_col = [x for x in flow_col if 'EpCAM' not in x]
# Spleen graphs
df = ad_mye.copy()
df = df.loc['spleen']
df[flow_col] = df[flow_col]*100
boxplot_titlestr = ' in Spleen'
pairplot_titlestr = 'Myeloid in Spleen'

create_plots(boxplot_titlestr, pairplot_titlestr, df, savedir, norm, dotplothue,
             mrti_col, flow_col, bxplt=True, pairplt=True)
#%%

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

#%% Compare thermal ablation amounts

df = ad_other.copy()
df = df.loc['tumor']
df = df.loc[(df['group2'] == 'fus') | (df['group2'] == 'combination')]
df['Hyper:TD240'] = df['Hyper_BV']/df['TD240_BV']

px = 1 / plt.rcParams['figure.dpi']
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(350 * px, 400 * px), dpi=100)

#fig = plt.figure(figsize=(9,8))
dfm = pd.melt(df, id_vars=['subject','group3'],value_vars=['TD240_BV','Hyper_BV'])
sns.boxplot(x='variable', y='value',data=dfm,hue='group3',palette=['#8E44AD', '#E67E22'])
plt.gca().set_ylabel('')
plt.gca().set_xlabel('')
plt.legend(bbox_to_anchor=(1.12, 0.55), loc='upper left', borderaxespad=0)
plt.show()


for col in mrti_col:
    #fig, ax = plt.subplot(1, 1, 1)
    sns.pointplot('group3', col, data=df, dodge=True, join=False, ci='sd', errwidth=.75,
                  capsize=.25, markers='s', color='black')
    sns.stripplot('group3', col, data=df, jitter=True, size=7, color='tab:blue')

    plt.show()


#%%

SOM = pd.read_csv('/Users/sarajohnson/Box Sync/Flow Cytometry Data/Johnson_Welm-Spectre/Spectre Analysis/Tumor/Tumor Live Results/rstudio-export/clustered.data.tsneDR.csv')

SOM = pd.read_csv('/Users/sarajohnson/Box Sync/Flow Cytometry Data/Johnson_Welm-Spectre/Spectre Analysis/Tumor/Tumor Live Results/annotation analysis-export/Annotated.data.csv')
