"""Export the trained glmnet-MLP weights to CSV files.

The CSV files (one per layer tensor, plus model_structure.csv) are the form in
which the trained network is embedded in the R package glmnetconf.

Input : results/eig{k}/best_model_glmnet_valid_eig{k}.pth
Output: results/eig{k}/model_params_glmnet_eig{k}/
Run from this directory:  python export_weights.py
"""
import os

import pandas as pd
import torch

EIG_K = 10   # final model (top/bottom 10 eigenvalues)

folder   = os.path.join("results", f"eig{EIG_K}")
ckpt     = torch.load(os.path.join(folder, f"best_model_glmnet_valid_eig{EIG_K}.pth"),
                      map_location="cpu")
out_dir  = os.path.join(folder, f"model_params_glmnet_eig{EIG_K}")
os.makedirs(out_dir, exist_ok=True)

for name, tensor in ckpt["model_state_dict"].items():
    pd.DataFrame(tensor.numpy()).to_csv(os.path.join(out_dir, f"{name}.csv"), index=False)

pd.DataFrame({
    "input_size"  : [ckpt["input_size"]],
    "hidden_sizes": [str(ckpt["hidden_sizes"])],
    "output_size" : [ckpt["output_size"]],
    "eig_k"       : [ckpt["eig_k"]],
}).to_csv(os.path.join(out_dir, "model_structure.csv"), index=False)
