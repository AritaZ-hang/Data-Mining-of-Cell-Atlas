# based on COLAB


from google.colab import drive
drive.mount('/content/drive')

import sys

# If True, will install via pypi, else will install from source
stable = True
IN_COLAB = "google.colab" in sys.modules

if IN_COLAB and stable:
    !pip install --quiet scvi-tools==0.14.6
elif IN_COLAB and not stable:
    !pip install --quiet --upgrade jsonschema
    !pip install --quiet git+https://github.com/yoseflab/scvi-tools@master#egg=scvi-tools[tutorials]

!pip install scanpy
!pip install bbknn
!pip install torch==1.10.1+cu102 torchvision==0.11.2+cu102 torchaudio===0.10.1+cu102 -f https://download.pytorch.org/whl/cu102/torch_stable.html
!pip install leidenalg
!pip install scikit-misc

import scvi
import anndata
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

adata_rna = sc.read('Yourpath/adata_rna.h5ad')
adata_atac_act = sc.read('Yourpath/adata_atac_act.h5ad')
adata_atac_act.obs.insert(adata_atac_act.obs.shape[1], 'Cell.type','Unknown')
adata_both = adata_rna.concatenate(adata_atac_act)

adata_both.obs

sc.pp.highly_variable_genes(adata_both, n_top_genes=4000, flavor = 'seurat_v3', batch_key = 'batch', subset = True)
scvi.data.setup_anndata(adata_both, labels_key='Cell.type', batch_key='batch')
model = scvi.model.SCVI(adata_both, gene_likelihood='nb', dispersion='gene-batch')
model.train()
model.save('Yourpath/SCVI_model', overwrite = True)
model = scvi.model.SCVI.load('Yourpath/SCVI_model', adata = adata_both)
lvae = scvi.model.SCANVI.from_scvi_model(model, unlabeled_category='Unknown', adata = adata_both)
lvae.train(max_epochs=100, n_samples_per_label=200)
lvae.save('Yourpath/lvae_model', overwrite = True)

adata_both.obs.insert(adata_both.obs.shape[1], 'predicted.labels', lvae.predict())
df = adata_both.obs
df = df[df.batch == '1']

df.to_csv('Yourpath/scvi_annotation.csv')
