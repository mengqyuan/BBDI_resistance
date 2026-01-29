########python###########
import sys, os, pickle
from collections import Counter
import importlib
from ipywidgets import widgets
import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
scriptpath = 'CaDRReS-Sc'
sys.path.append(os.path.abspath(scriptpath))
from cadrres_sc import pp, model, evaluation, utility
model_dir = 'CaDRReS-Sc/example_result/'
obj_function = widgets.Dropdown(options=['cadrres-wo-sample-bias', 'cadrres-wo-sample-bias-weight'], description='Objetice function')
model_spec_name = obj_function.value
model_file = model_dir + '{}_param_dict.pickle'.format(model_spec_name)
cadrres_model = model.load_model(model_file)
gdsc_gene_exp_df = pd.read_csv('/data/home/yuanmengqin/CaDRReS-Sc/data/GDSC/GDSC_exp.tsv', sep='\t', index_col=0)
gdsc_gene_exp_df = gdsc_gene_exp_df.groupby(gdsc_gene_exp_df.index).mean()
gdsc_norm_exp_df, gdsc_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gdsc_gene_exp_df)
ess_gene_list = utility.get_gene_list('CaDRReS-Sc/data/essential_genes.txt')

gene_exp_df = pd.read_csv('SUM159_exp.tsv', sep='\t', index_col=0)

cell_line_log2_mean_fc_exp_df, cell_line_mean_exp_df = pp.gexp.normalize_log2_mean_fc(gene_exp_df)
ess_gene_list = utility.get_gene_list('CaDRReS-Sc/data/essential_genes.txt')

test_kernel_df = pp.gexp.calculate_kernel_feature(cell_line_log2_mean_fc_exp_df, gdsc_norm_exp_df, ess_gene_list)
#Calculating kernel features based on 469 common genes
print("Dataframe shape:", test_kernel_df.shape, "\n")
test_kernel_df.head(2)
print('Predicting drug response using CaDRReS: {}'.format(model_spec_name))
pred_df, P_test_df= model.predict_from_model(cadrres_model, test_kernel_df, model_spec_name)
print('done!')
print('Saving ' + model_dir + '{}_test_pred.csv'.format(model_spec_name))
pred_df.to_csv("log2_IC50.pred.csv") #Predicted drug response (log2 IC50)