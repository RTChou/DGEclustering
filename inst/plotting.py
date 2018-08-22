#!/usr/bin/env python3

import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
from scipy.stats.mstats import mquantiles
from scipy.stats import beta

def scatter_plot(file_path_1, file_path_2, plot_out_dir='./', dat_out_dir='./', x_threshold=0.05, y_threshold=0.05, adj_pvalue=True, for_cluster_plot=False, return_sig_plot=False):
    
    plt.close()
    dataset = pd.read_table(file_path_1)
    dataset_2 = pd.read_table(file_path_2)
    merged_set = dataset.merge(dataset_2, left_on=dataset.columns[0], right_on=dataset_2.columns[0])
    
    if adj_pvalue == True:
        sig_vs_sig = merged_set[(merged_set['padj_x'] < x_threshold) & (merged_set['padj_y'] < y_threshold)]
        sig_vs_NS = merged_set[(merged_set['padj_x'] < x_threshold) & (merged_set['padj_y'] >= y_threshold)]
        NS_vs_sig = merged_set[(merged_set['padj_x'] >= x_threshold) & (merged_set['padj_y'] < y_threshold)]
        NS_vs_NS = merged_set[(merged_set['padj_x'] >= x_threshold) & (merged_set['padj_y'] >= y_threshold)]
        non_NA_set = merged_set[(merged_set.isna()['padj_x'] == False) & (merged_set.isna()['padj_y'] == False)]
    
    else:
        sig_vs_sig = merged_set[(merged_set['pvalue_x'] < x_threshold) & (merged_set['pvalue_y'] < y_threshold)]
        sig_vs_NS = merged_set[(merged_set['pvalue_x'] < x_threshold) & (merged_set['pvalue_y'] >= y_threshold)]
        NS_vs_sig = merged_set[(merged_set['pvalue_x'] >= x_threshold) & (merged_set['pvalue_y'] < y_threshold)]
        NS_vs_NS = merged_set[(merged_set['pvalue_x'] >= x_threshold) & (merged_set['pvalue_y'] >= y_threshold)]
        non_NA_set = merged_set[(merged_set.isna()['pvalue_x'] == False) & (merged_set.isna()['pvalue_y'] == False)]  

    sig_discordant = sig_vs_sig[((sig_vs_sig['log2FoldChange_x'] < 0) & (sig_vs_sig['log2FoldChange_y'] > 0)) | 
             ((sig_vs_sig['log2FoldChange_x'] > 0) & (sig_vs_sig['log2FoldChange_y'] < 0))]    
    sig_concordant = sig_vs_sig[((sig_vs_sig['log2FoldChange_x'] >= 0) & (sig_vs_sig['log2FoldChange_y'] >=0)) |
             ((sig_vs_sig['log2FoldChange_x'] <= 0) & (sig_vs_sig['log2FoldChange_y'] <= 0))]

    fig = plt.figure(figsize=(18, 18))
    ax = fig.add_subplot(111)
    if for_cluster_plot == False and return_sig_plot == False:
        g4 = ax.scatter(NS_vs_NS['log2FoldChange_x'], NS_vs_NS['log2FoldChange_y'], s=9, c='grey', alpha=0.3)
        g3 = ax.scatter(NS_vs_sig['log2FoldChange_x'], NS_vs_sig['log2FoldChange_y'], s=9, c=(31 / 255., 119 / 255., 180 / 255.), alpha=0.6)
        g2 = ax.scatter(sig_vs_NS['log2FoldChange_x'], sig_vs_NS['log2FoldChange_y'], s=9, c=(255 / 255., 127 / 255., 14 / 255.), alpha=0.6)
        g1 = ax.scatter(sig_vs_sig['log2FoldChange_x'], sig_vs_sig['log2FoldChange_y'], s=15, c=(214 / 255., 39 / 255., 40 / 255.), alpha=1.0)

    filename_1 = re.search(r".+\/(.+).tsv", file_path_1).group(1)
    filename_2 = re.search(r".+\/(.+).tsv", file_path_2).group(1)
    xtitle = filename_1.replace('_', ' ')
    xtitle = xtitle.replace('.', ' ')
    ytitle = filename_2.replace('_', ' ')
    ytitle = ytitle.replace('.', ' ')
    ax.set_xlim(min(non_NA_set['log2FoldChange_x'].min(), non_NA_set['log2FoldChange_y'].min()) - 0.5, 
            max(non_NA_set['log2FoldChange_x'].max(), non_NA_set['log2FoldChange_y'].max()) + 0.5)
    ax.set_ylim(min(non_NA_set['log2FoldChange_x'].min(), non_NA_set['log2FoldChange_y'].min()) - 0.5, 
            max(non_NA_set['log2FoldChange_x'].max(), non_NA_set['log2FoldChange_y'].max()) + 0.5)
    ax.axvline(x=0, linestyle='dotted', color='grey')
    ax.axhline(y=0, linestyle='dotted', color='grey')
    ax.set_title('(' + xtitle + ') vs (' + ytitle + ') (gene number=' + str(merged_set.shape[0]) + ')', fontweight='bold', fontsize=16, y=1.02)
    ax.set_xlabel(xtitle + u' log\u2082 fold change', fontsize=15)
    ax.set_ylabel(ytitle + u' log\u2082 fold change', fontsize=15)
    
    if for_cluster_plot == False and return_sig_plot == False:
        ax.legend((g1, g2, g3, g4),(
            'sig vs sig (' + str(sig_vs_sig.shape[0]) + ')',
            'sig vs NS (' + str(sig_vs_NS.shape[0]) + ')',
            'NS vs sig (' + str(NS_vs_sig.shape[0]) + ')',
            'NS vs NS (' + str(NS_vs_NS.shape[0]) + ')'),
            markerscale=1)
    anchored_text = AnchoredText('# of sig vs sig in II and IV: ' + str(sig_discordant.shape[0]), loc=3)
    anchored_text.patch.set(color='red', alpha=0.3)
    ax.add_artist(anchored_text)    

    if for_cluster_plot == False and return_sig_plot == False:
        fig.savefig(plot_out_dir + '/' + filename_1 + '_vs_' + filename_2 + '_scatter_plot.png')
        if sig_discordant.shape[0] > 0:
            discordant_path = dat_out_dir + '/' + filename_1 + '_vs_' + filename_2 + '_disagreeing_genes.tsv'
            sig_discordant.to_csv(discordant_path, sep='\t', index=False)
        else:
            discordant_path = None

        if sig_concordant.shape[0] > 0:
            concordant_path = dat_out_dir + '/' + filename_1 + '_vs_' + filename_2 + '_agreeing_genes.tsv'
            sig_concordant.to_csv(concordant_path, sep='\t', index=False)
        else:
            concordant_path = None

        return {'discordant_path': discordant_path, 'concordant_path': concordant_path}
    
    elif for_cluster_plot == True:
        return {'plot': ax, 'discordant': sig_discordant, 'concordant': sig_concordant}   
    
    else:
        g4 = ax.scatter(NS_vs_NS['log2FoldChange_x'], NS_vs_NS['log2FoldChange_y'], s=9, c='grey', alpha=0.3)
        g3 = ax.scatter(NS_vs_sig['log2FoldChange_x'], NS_vs_sig['log2FoldChange_y'], s=9, c='grey', alpha=0.6)
        g2 = ax.scatter(sig_vs_NS['log2FoldChange_x'], sig_vs_NS['log2FoldChange_y'], s=9, c='grey', alpha=0.6)
        g1 = ax.scatter(sig_vs_sig['log2FoldChange_x'], sig_vs_sig['log2FoldChange_y'], s=15, c=(214 / 255., 39 / 255., 40 / 255.), alpha=1.0)
        ax.legend((g1,), ('sig vs sig (' + str(sig_vs_sig.shape[0]) + ')',), markerscale=1)
        return ax


def fish_plot(file_path_1, file_path_2, output_dir):
  
    plt.close()  
    dataset = pd.read_table(file_path_1)
    dataset_2 = pd.read_table(file_path_2)
    merged_set = dataset.merge(dataset_2, left_on=dataset.columns[0], right_on=dataset_2.columns[0])
    merged_set['-log10_pvalue_x'] = - np.sign(merged_set['log2FoldChange_x']) * np.sign(merged_set['log2FoldChange_y']) * np.log10(merged_set['pvalue_x'])
    merged_set['-log10_pvalue_y'] = - np.sign(merged_set['log2FoldChange_x']) * np.sign(merged_set['log2FoldChange_y']) * np.log10(merged_set['pvalue_y'])

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.scatter(merged_set['-log10_pvalue_x'], merged_set['-log10_pvalue_y'], s=3, c=(31 / 255., 119 / 255., 180 / 255.), alpha=0.5)

    filename_1 = re.search(r".+\/(.+).tsv", file_path_1).group(1)
    filename_2 = re.search(r".+\/(.+).tsv", file_path_2).group(1)
    xtitle = filename_1.replace('_', ' ')
    xtitle = xtitle.replace('.', ' ')
    ytitle = filename_2.replace('_', ' ')
    ytitle = ytitle.replace('.', ' ')
    #ax.set_xlim(min(merged_set['-log10_pvalue_x'].min(), merged_set['-log10_pvalue_y'].min()) - 0.5,
    #        max(merged_set['-log10_pvalue_x'].max(), merged_set['-log10_pvalue_y'].max()) + 0.5)
    #ax.set_ylim(min(merged_set['-log10_pvalue_x'].min(), merged_set['-log10_pvalue_y'].min()) - 0.5,
    #        max(merged_set['-log10_pvalue_x'].max(), merged_set['-log10_pvalue_y'].max()) + 0.5)
    ax.set_xlim(-40, 40)
    ax.set_ylim(-40, 40)
    ax.axvline(x=0, linestyle='dotted', color='grey')
    ax.axhline(y=0, linestyle='dotted', color='grey')
    ax.set_title('(' + xtitle + ') vs (' + ytitle + ') (gene number=' + str(merged_set.shape[0]) + ')', fontweight='bold', fontsize=16, y=1.02)
    ax.set_xlabel(xtitle + u' -log\u2081\u2080 pvalue', fontsize=15)
    ax.set_ylabel(ytitle + u' -log\u2081\u2080 pvalue', fontsize=15)

    fig.savefig(output_dir + '/' + filename_1 + '_vs_' + filename_2 + '_fish_plot.png')


def qq_plot(output_dir, file_path=None, dataset=None):
    
    plt.close()
    if dataset is None:
        dataset = pd.read_table(file_path)
    dataset = dataset[dataset.isna()['pvalue'] == False]
    dataset['-log10_pvalue'] = - np.log10(dataset['pvalue'])
    gene_size = dataset.shape[0]

    exp = np.concatenate([np.arange(100) / gene_size, np.logspace(-np.log10(gene_size) + 2, 0, 200)])
    obs = mquantiles(dataset['pvalue'], prob=exp, alphap=0, betap=1)

    lower = list()
    upper = list()
    for i in range(0, len(exp)):
        CI_values = beta.interval(0.95, gene_size * exp[i], gene_size - gene_size * exp[i])
        lower.append(CI_values[0])
        upper.append(CI_values[1])

    exp = -np.log10(exp)
    obs = -np.log10(obs)
    up = -np.log10(lower)
    low = -np.log10(upper)

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax.fill_between(exp, up, low, color='grey', alpha=0.5)
    ax.set_xlim(np.nanmin(exp[exp != -np.inf]), np.nanmax(exp[exp != np.inf]) + 0.1)
    ax.set_ylim(np.nanmin(obs[obs != -np.inf]), max(np.nanmax(obs[obs != np.inf]), np.nanmax(up[up != np.inf])) + 0.5)
    ax.plot(ax.get_xlim(), ax.get_xlim(), linestyle='--', color='black')
    ax.scatter(exp, obs, s=3, c=(31 / 255., 119 / 255., 180 / 255.))

    if file_path is not None:
        filename = re.search(r".+\/(.+).tsv", file_path).group(1).replace('_', ' ')
    else:
        filename = 'Null'
    ax.set_title(filename + ' QQ-Plot', fontweight='bold', fontsize=16, y=1.02)
    ax.set_xlabel('expected -log\u2081\u2080 pvalue', fontsize=15)
    ax.set_ylabel('observed -log\u2081\u2080 pvalue', fontsize=15)
    #ax.axvline(x=0, linestyle='dotted', color='grey')
    #ax.axhline(y=0, linestyle='dotted', color='grey')

    fig.savefig(output_dir + '/' + filename + '_qq_plot.png')

