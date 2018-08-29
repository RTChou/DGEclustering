#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import plotting
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.patches import Ellipse
from scipy.stats import linregress
import math
import warnings

def main():
    parser = argparse.ArgumentParser(description="This script is for clustering visulization of RNA-seq paired DGE files")
    parser.add_argument('-d', '--out_dir', required=True, help='output directory for resulting plots')
    parser.add_argument('-l', '--list', nargs='+', required=True, help='a list of filepaths. e.g. ./multidimension.py -l path1 path2 path3')
    parser.add_argument('-n1', 'x_file_number', default=0, type=int, help='file number in list for x axis. index starts at 0')
    parser.add_argument('-n2', 'y_file_number', default=0, type=int, help='file number in list for y axis. index starts at 0')
    parser.add_argument('-g', '--gene_col', required=True, type=str, help='gene ID column name') 
    parser.add_argument('-r', '--clustering_result', required=True, help='file path of the clustering result')
    parser.add_argument('-m', '--MCA_result', help='file path of the MCA result')
  #  parser.add_argument('-x', '--x_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot x axis')
  #  parser.add_argument('-y', '--y_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot y axis')
  #  parser.add_argument('-a', '--adj_pvalue', default=1, type=int, help='whether to use adjusted pvalue or pvalue. 1 as True, 0 as False')
  #  parser.add_argument('-s', '--sig_data', default='all', help='one of \'dis\', \'con\', \'all\' (dis + con), or \'multi\'')
    parser.add_argument('-c', '--color', default='brg', help='cmap color for visualization')
    args = parser.parse_args()
    
    warnings.filterwarnings('ignore') # ignore runtime warnings
    
    # generate sig data-only plot
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    ax = plotting.scatter_plot(file_paths=args.list, x_file_number=args.x_file_number, y_file_number=args.y_file_number, x_threshold=args.x_threshold, y_threshold=args.y_threshold, adj_pvalue=args.adj_pvalue, return_sig_plot=True)
    plt.savefig(args.out_dir + '/sig_data.png')

    # prepare for cluster plotting
    plt.close()
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)
    temp = plotting.scatter_plot(file_paths=args.list, x_file_number=args.x_file_number, y_file_number=args.y_file_number, x_threshold=args.x_threshold, y_threshold=args.y_threshold, adj_pvalue=args.adj_pvalue, for_cluster_plot=True)
    ax = temp['plot']
    
    clustered_dat = pd.read_table(args.clustering_result)   
    ## merge with file1 file2
    sig = sig.merge(clustered_dat, left_on=('0_' + args.gene_col), right_on=clustered_dat.columns[0])
 
    # plot clusters from significant data
    i = 1
    group_list = list()
    while ('Group.' + str(i)) in sig['ind'].tolist():
        group_list.append(sig[sig['ind'] == 'Group.' + str(i)].reset_index(drop=True))
        i += 1   
    cmap = cm.get_cmap(args.color)(np.linspace(0, 1.0, len(group_list)))
    if args.MCA_result is not None:
        MCA = pd.read_table(args.MCA_result)
    log2FoldChange_x = str(x_file_number) + '_log2FoldChange'
    log2FoldChange_y = str(y_file_number) + '_log2FoldChange'

    # main plot
    groups = list()
    group_names = list()  
    for index, group in enumerate(group_list):
        color = cmap[index]
        x = group[log2FoldChange_x]
        y = group[log2FoldChange_y]

        # plot the centroid for each cluster
        cx = np.mean(x)
        cy = np.mean(y)
        ax.scatter(cx, cy, s=9, facecolors='none', edgecolors='black', alpha=0.5)
        ax.annotate(str(index + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01))

        # plot data points
        groups.append(ax.scatter(x, y, s=10, c=([color,] * group.shape[0]), alpha=0.6))
        group_names.append('cluster ' + str(index + 1) + ' (' + str(group.shape[0]) + ')')

        # plot lines from centroid to points
        for row in range(0, group.shape[0]):
            ax.plot([cx, x[row]], [cy, y[row]], 
                    linestyle='dashed', linewidth=0.5, c=color, alpha=0.3)
        
        # encircle the group points
        slope, intercept, r_value, p_value, std_err = linregress(x,y)
        angle = math.degrees(math.atan(slope))
        if np.std(x) < np.std(y):
                    angle = angle + 90
        ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5, angle=angle, ec='k', fc=color, alpha=0.2)
        ax.add_patch(ell)
    
    ax.legend(groups, group_names, markerscale=1)
    # plt.show()
    plt.savefig(args.out_dir + '/cluster_all.png')

    # subplots
    for j in range(0, len(group_list)):
        plt.close()
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111)
        temp = plotting.scatter_plot(file_paths=args.list, x_file_number=args.x_file_number, y_file_number=args.y_file_number, x_threshold=args.x_threshold, y_threshold=args.y_threshold, adj_pvalue=args.adj_pvalue, for_cluster_plot=True)
        ax = temp['plot']

        for index, group in enumerate(group_list):
            color = cmap[index]
            x = group[log2FoldChange_x]
            y = group[log2FoldChange_y]

            if index != j:
                # plot the centroid for each cluster
                cx = np.mean(x)
                cy = np.mean(y)
                ax.scatter(cx, cy, s=3, facecolors='none', edgecolors='grey', alpha=0.3)
                ax.annotate(str(index + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01), color='grey')

                # plot lines from centroid to points
                for row in range(0, group.shape[0]):
                    ax.plot([cx, x[row]], [cy, y[row]],
                            linestyle='dashed', linewidth=0.5, c='grey', alpha=0.3)
                 
                # plot data points
                ax.scatter(x, y, s=5, c=('grey'), alpha=0.3)

                # encircle the group points
                slope, intercept, r_value, p_value, std_err = linregress(x,y)
                angle = math.degrees(math.atan(slope))
                if np.std(x) < np.std(y):
                    angle = angle + 90
                ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5,
                        angle=angle, ec='#7f7f7f', fc='grey', alpha=0.1)
                ax.add_patch(ell)

        
        # plot the centroid for cluster j
        color = cmap[j]
        x = group_list[j][log2FoldChange_x]
        y = group_list[j][log2FoldChange_y]
        
        cx = np.mean(x)
        cy = np.mean(y)
        ax.scatter(cx, cy, s=9, facecolors='none', edgecolors='grey')
        ax.annotate(str(j + 1), xy=(cx, cy), xytext=(cx + 0.01, cy + 0.01))
                
        # plot lines from centroid to cluster j points
        for row in range(0, group_list[j].shape[0]):
            ax.plot([cx, x[row]], [cy, y[row]], linestyle='dashed', linewidth=0.5, c=color)
        
        # plot data points
        g1 = ax.scatter(x, y, s=20, c=([color,] * group.shape[0]), edgecolors='black', linewidth=0.8)
        g1_name = 'cluster ' + str(j + 1) + ' (' + str(group_list[j].shape[0]) + ')'
                 
        # encircle the cluster j group points
        slope, intercept, r_value, p_value, std_err = linregress(x,y)
        angle = math.degrees(math.atan(slope))
        if np.std(x) < np.std(y):
             angle = angle + 90 
        ell = Ellipse(xy=[cx, cy], width=np.std(x) * 5, height=np.std(y) * 5, angle=angle, ec='k', fc=color, alpha=0.2)
        ax.add_patch(ell)

        ax.legend((g1,), (g1_name,), markerscale=1)
        plt.savefig(args.out_dir + '/cluster_' + str(j + 1) + '.png')


if __name__ == '__main__':
    main()

