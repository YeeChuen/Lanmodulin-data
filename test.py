import os
import numpy as np
import matplotlib.pyplot as plt

def heatmapPlot(data, save_name):
    plt.figure(figsize=(10, 8))
    plt.imshow(data, cmap='viridis', origin='lower', interpolation='nearest')
    plt.colorbar(label='RMSD')
    plt.title('RMSD Heatmap')
    plt.xlabel('PDB File Index')
    plt.ylabel('PDB File Index')
    plt.savefig(save_name)


data0 = [[0,1,2], [1,0,3], [2,3,0]]
data1 = [[0,5,2], [5,0,3], [2,3,0]]
data2 = [[0,1,4], [1,0,3], [4,3,0]]

for i, data in enumerate([data0, data1, data2]):
    heatmapPlot(data, f'heatmap_{i}')

