# Author: Yee Chuen Teoh
# python build_barchart.py

#______
# imports

#plot library
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import seaborn as sns
import numpy as np
from matplotlib_venn import venn2

#______
# functions
# check for unique seq
def checkUnique():
    motif_dict = {}
    seq_occurence = {}
    with open('EFhands.fasta', 'r') as f:
        contents = f.readlines()
        contents = [x.replace("\n", "") for x in contents]

        for i in range(0, len(contents), 2):
            motif = contents[i].split("|")[-1]
            seq = contents[i + 1]
            if motif not in motif_dict: motif_dict[motif] = []
            if seq not in seq_occurence: seq_occurence[seq] = 0
            motif_dict[motif].append(seq)
            seq_occurence[seq] += 1
    
    #print(motif_dict) # <-- development and debug purposes
    #print(seq_occurence) # <-- development and debug purposes

    total = sum([seq_occurence[key] for key in seq_occurence])
    #print(total) # <-- development and debug checker
    # 2942 because we have 623 lanmodulin, each lanmodulin has 4 hands.
    if total != 2492: raise ValueError(f"\nTotal EF hands does not match 623 * 4 = 2492,\nData has {total} EF")
        
    data = [len(motif_dict[key]) for key in motif_dict]

    barChartPlot(data = data, labels = list(motif_dict.keys()), 
                y_label = "EF hands", title = "Number of EF hands to motifs", 
                save_name = "bar_chart_ef_to_motif", y_limit = None, bar_color = "#0675b5")

    x_order = list(motif_dict.keys())
    #print(x_order)
    #print(motif_dict)
    #print(seq_occurence)

    # we dont overwrite data so we know how many slot we need for x axis
    data_stacked = []
    for seq_key in seq_occurence: # <-- unique sequence
        data_row = [0 for _ in range(len(x_order))]
        for motif_key in motif_dict:
            if seq_key in motif_dict[motif_key]:
                data_row[x_order.index(motif_key)] = seq_occurence[seq_key]
        data_stacked.append(data_row)

    stackedBarChartPlot(data = data_stacked, labels = x_order, 
                y_label = "EF hands", title = "Number of Unique EF hands in motifs", 
                save_name = "bar_chart_unique_ef_to_motif", y_limit = None)

'''
Stacked plot bar chart.
Modified from VE2_plot_helper.py 
'''
def stackedBarChartPlot(data = [[5, 6, 3], [2, 1, 3], [4, 2, 2], [2, 4, 9]], labels = ['x1', 'x2', 'x3'], 
                y_label = "y_label", title = "Bar chart title", 
                save_name = "test_barChart", y_limit = None):
    fig, ax = plt.subplots()
    
    y_pos = np.arange(len(labels))
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False

    # stacked here.
    prev = [0 for _ in range(len(data[0]))]
    for d in data:
        plt.bar(labels, d, bottom=prev)
        prev = [prev[i] + d[i] for i in range(len(prev))]
        #print(d)

    plt.xticks(y_pos, labels, rotation = 90)
    plt.ylabel(y_label)
    plt.title(title)

    if y_limit:
        plt.ylim([y_limit[0], y_limit[1] + math.floor(y_limit[1] * 0.1)])
        
    fig.tight_layout()
    plt.savefig(save_name)
    plt.close()
    pass

'''
Basic plot bar chart.
Taken from VE2_plot_helper.py 
from path: /work/ratul1/chuen/viral_escape/pyRosetta_script
'''
def barChartPlot(data = [5, 6, 3], labels = ['x1', 'x2', 'x3'], y_label = "y_label", title = "Bar chart title", 
                save_name = "test_barChart", y_limit = None, bar_color = "#0675b5"):
    fig, ax = plt.subplots()
    
    y_pos = np.arange(len(labels))
    
    plt.rcParams['axes.spines.right'] = False
    plt.rcParams['axes.spines.top'] = False

    color = [bar_color for _ in range(len(labels))]
    plt.bar(y_pos, data, align='center', alpha=0.5, color = color)
    plt.xticks(y_pos, labels, rotation = 90)
    plt.ylabel(y_label)
    plt.title(title)

    if y_limit:
        plt.ylim([y_limit[0], y_limit[1] + math.floor(y_limit[1] * 0.1)])
        
    fig.tight_layout()
    plt.savefig(save_name)
    plt.close()
    pass
    
#______
# main
def main():
    checkUnique()
    pass

if __name__ == "__main__":
    main()
