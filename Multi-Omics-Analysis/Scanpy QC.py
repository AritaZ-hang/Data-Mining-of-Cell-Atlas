import scanpy as sc
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re

working_dir = r'Yourpath'
os.chdir(working_dir)

if not os.path.exists(working_dir + '/temp_data'):
    os.mkdir(working_dir + '/temp_data')
if not os.path.exists(working_dir + '/dealt_data'):
    os.mkdir(working_dir + '/dealt_data')

temp_path = working_dir + '/temp_data'
save_path = working_dir + '/dealt_data'

def MatrixReverse(file_name):
    origin_matrix = pd.read_table(file_name, sep = ' ')
    reverse_matrix = pd.DataFrame(origin_matrix.values.T, index = origin_matrix.columns, columns = origin_matrix.index)
    pattern = re.compile(r'^.*?_')
    save_name = pattern.match(file_name).group(0)
    reverse_matrix.to_csv(temp_path + '/' + save_name + 'matrix.csv')

def preprocessing(csv_name):
    adata = sc.read_csv(csv_name)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_counts = 500)
    sc.pp.filter_genes(adata, min_cells = 3)
    adata.var['mt'] = adata.var_names.str.startswith('mt')
    sc.pp.calculate_qc_metrics(adata, qc_vars = ['mt'], percent_top=None, log1p=None, inplace=True)
    adata = adata[adata.obs.n_genes_by_counts < 2500,:]
    adata = adata[adata.obs.pct_counts_mt < 5,:]
    adata = adata[:,~adata.var_names.str.startswith('mt-')]
    
    pattern = re.compile(r'^.*?_')
    save_name = pattern.match(csv_name).group(0)
    adata.write(temp_path + '/' + save_name + 'matrix.h5ad')
