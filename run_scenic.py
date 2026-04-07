import os
import numpy as np
import scanpy as sc
import loompy as lp
import subprocess

print(os.getcwd())
print(os.listdir(os.getcwd()))

x = sc.read_csv("for.scenic.data.csv")

row_attrs = {"Gene": np.array(x.var_names)}
col_attrs = {"CellID": np.array(x.obs_names)}

lp.create("sample.loom", x.X.transpose(), row_attrs, col_attrs)

print("LOOM file created: sample.loom")

dir_path = "/home/lipeixin27/data/index_genome/cisTarget_databases/"

tfs = os.path.join(dir_path, "hs_hgnc_tfs.txt")
feather = os.path.join(dir_path, "hg19-tss-centered-10kb-10species.mc9nr.genes_vs_motifs.rankings.feather")
tbl = os.path.join(dir_path, "motifs-v9-nr.hgnc-m0.001-o0.0.tbl")

input_loom = "./sample.loom"

for f in [tfs, feather, tbl]:
    if not os.path.exists(f):
        raise FileNotFoundError(f"{f} not found")

print("All reference files exist.")

# GRN
subprocess.run([
    "pyscenic", "grn",
    "--num_workers", "10",
    "--output", "adj.sample.tsv",
    "--method", "grnboost2",
    input_loom,
    tfs
], check=True)

# CTX
subprocess.run([
    "pyscenic", "ctx",
    "adj.sample.tsv",
    feather,
    "--annotations_fname", tbl,
    "--expression_mtx_fname", input_loom,
    "--mode", "dask_multiprocessing",
    "--output", "reg.csv",
    "--num_workers", "20",
    "--mask_dropouts"
], check=True)

# AUCELL
subprocess.run([
    "pyscenic", "aucell",
    input_loom,
    "reg.csv",
    "--output", "out_SCENIC.loom",
    "--num_workers", "10"
], check=True)

print("SCENIC pipeline finished!")
