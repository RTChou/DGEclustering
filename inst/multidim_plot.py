#!/usr/bin/env python3

import argparse
import warnings

def main():
    parser = argparse.ArgumentParser(description="This script is for plotting and generating datasets for multi-dimensional clustering")
    parser.add_argument('-f', '--file_paths', nargs='+', required=True, help='a list of filepaths. e.g. ./multidimension.py -l path1 path2 path3')
    parser.add_argument('-n1', 'x_file_number', default=0, type=int, help='file number in list for x axis. index starts at 0')
    parser.add_argument('-n2', 'y_file_number', default=1, type=int, help='file number in list for y axis. index starts at 0')
    parser.add_argument('-p', '--plot_dir', required=True, help='plotting directory')
    parser.add_argument('-d', '--dat_dir', required=True, help='significant data directory')
    parser.add_argument('-o', '--output', help='output file name')
    parser.add_argument('-x', '--x_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot x axis')
    parser.add_argument('-y', '--y_threshold', default=0.05, type=float, help='(adjusted) pvalue for scatter plot y axis')
    parser.add_argument('-a', '--adj_pvalue', default=1, type=int, help='whether to use adjusted pvalue or pvalue. 1 as True, 0 as False')
    args = parser.parse_args()

    warnings.filterwarnings('ignore') # ignore runtime warnings
    temp = plotting.scatter_plot(file_paths=args.file_paths, x_file_number=args.x_file_number, y_file_number=args.y_file_number, plot_out_dir=args.plot_dir, dat_out_dir=args.dat_dir, x_threshold=args.x_threshold, y_threshold=args.y_threshold, adj_pvalue=args.adj_pvalue, return_sig_plot=True, out_file_name=args.output)
    
if __name__ == '__main__':
    main()

