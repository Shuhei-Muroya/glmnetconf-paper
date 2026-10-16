"""Section 3.2.2 / Appendix C: training of glmnet-MLP.

- Bayesian optimisation (Optuna + BoTorchSampler, EI, 500 trials x 100 epochs)
  over the number of hidden layers {1,2,3}, units {1..64} and learning rate.
- Final training: up to 1000 epochs, early stopping (patience 100),
  minibatch 256, Swish activation, 80/10/10 train/valid/test split.
- Run for the top/bottom k eigenvalue features, k in {10, 5, 1}
  ("Number of eigenvalue features" in Appendix C). k = 10 is the final model
  (hidden layers 49-63-32), shipped in the R package glmnetconf.

Input : ../01_summary_dataset/summary_dataset.csv.gz
Output: results/eig{k}/  (model checkpoint, scalers, Optuna study, plots)
        results/eigenvalue_ablation_summary.csv, results/eigenvalue_scaled_rmse.csv
Run from this directory:  python train_glmnet_mlp.py
"""
# ========= Comparison of the number of eigenvalue features: top/bottom 10 -> 5 -> 1 =========
# Bayesian optimisation and epoch selection use VALID; the final evaluation uses TEST
# the acquisition function is defined in acquisition.py
from acquisition import EI_candidates_func

import optuna
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.metrics import mean_squared_error
from sklearn.preprocessing import StandardScaler
from torch.utils.data import DataLoader, TensorDataset
import pandas as pd
import os
import random
import copy
import matplotlib.pyplot as plt

# =========================
# Settings
# =========================
foldername_base = "results"
filename        = "../01_summary_dataset/summary_dataset"
batch_size      = 256
eps             = 0          # to avoid log(0) (set a positive value if needed)
n_bo_trials     = 500
max_epochs      = 1000
patience        = 100        # early stopping: stop if the VALID RMSE does not improve for this many epochs
seed            = 1

# numbers of eigenvalue features to run, in this order
# e.g. [10, 5, 1]: top/bottom 10, then top/bottom 5, then top/bottom 1
EIG_K_LIST = [10, 5, 1]

y_colname = ['time_glmnet', 'SPE_glmnet']


def make_x_colname(eig_k):
    """Input columns using the top/bottom eig_k eigenvalue features."""
    return (
        ['n', 'p']
        + [f'sample_ev_top{i}'    for i in range(1, eig_k + 1)]
        + [f'sample_ev_bottom{i}' for i in range(1, eig_k + 1)]
        + ['nlambda', 'thresh']
    )


# GPU
# gpu_id: index of the GPU (see nvidia-smi)
gpu_id = 1
if torch.cuda.is_available() and gpu_id < torch.cuda.device_count():
    device = torch.device(f'cuda:{gpu_id}')
else:
    device = torch.device('cpu')
print(f"Using device: {device}")
if device.type == 'cuda':
    print(f"  GPU {gpu_id}: {torch.cuda.get_device_name(gpu_id)}")

os.makedirs(foldername_base, exist_ok=True)

# =========================
# Model
# =========================
class Swish(nn.Module):
    def forward(self, x):
        return x * torch.sigmoid(x)


class Net(nn.Module):
    def __init__(self, input_size, hidden_sizes, output_size):
        super().__init__()
        layers = []
        in_size = input_size
        for h in hidden_sizes:
            layers += [nn.Linear(in_size, h), Swish()]
            in_size = h
        layers.append(nn.Linear(in_size, output_size))
        self.model = nn.Sequential(*layers)

    def forward(self, x):
        return self.model(x)


# =========================
# Fix the seed
# =========================
def set_seed(s):
    torch.manual_seed(s)
    torch.cuda.manual_seed_all(s)
    np.random.seed(s)
    random.seed(s)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# =========================
# Plot of predictions vs observations (returns the RMSE of each column)
# =========================
def plot_predictions_vs_observations(y_df, predictions, folder_name, name=""):
    os.makedirs(folder_name, exist_ok=True)
    rmse_dict = {}
    for i, col in enumerate(y_df.columns):
        plt.figure()
        plt.scatter(predictions[:, i], y_df[col], label='Data Points')
        min_val = min(np.min(predictions[:, i]), y_df[col].min())
        max_val = max(np.max(predictions[:, i]), y_df[col].max())
        plt.plot([min_val, max_val], [min_val, max_val], label='y=x')
        rmse = np.sqrt(mean_squared_error(y_df[col], predictions[:, i]))
        rmse_dict[col] = rmse
        plt.xlabel('Prediction')
        plt.ylabel('Observation')
        plt.title(f'Prediction vs Observation: {col}')
        plt.legend(loc='upper left')
        plt.text(max_val, min_val, f'RMSE: {rmse:.4f}', ha='right', va='bottom', fontsize=12)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.xlim(min_val, max_val)
        plt.ylim(min_val, max_val)
        safe_col = col.replace(".", "_").replace(" ", "_")
        suffix = f"_{name}" if name else ""
        plt.savefig(os.path.join(folder_name, f'plot_{safe_col}{suffix}.pdf'), bbox_inches='tight')
        plt.close()
    return rmse_dict


# =========================
# Read the data and split the rows (independent of eig_k, done once)
# =========================
csv_file_path = f"{filename}.csv.gz"   # pandas reads gzip directly
combined_df = pd.read_csv(csv_file_path)
n = combined_df.shape[0]

n_train = round(n * 0.8)
n_valid = round(n * 0.1)
n_test  = n - n_train - n_valid

train_df = combined_df.sample(n=n_train, random_state=1)
remaining_df = combined_df.drop(train_df.index)
valid_df = remaining_df.sample(n=n_valid, random_state=2)
test_df  = remaining_df.drop(valid_df.index)

print(f"Train: {train_df.shape}, Valid: {valid_df.shape}, Test: {test_df.shape}")

# y does not depend on eig_k, so it is processed once
y_train = train_df[y_colname]
y_valid = valid_df[y_colname]
y_test  = test_df[y_colname]

y_train_log = np.log(y_train + eps)
y_valid_log = np.log(y_valid + eps)
y_test_log  = np.log(y_test  + eps)

scaler_y = StandardScaler()
y_train_scaled = scaler_y.fit_transform(y_train_log)
y_valid_scaled = scaler_y.transform(y_valid_log)
y_test_scaled  = scaler_y.transform(y_test_log)

num_q = y_train.shape[1]

y_train_tensor = torch.tensor(y_train_scaled, dtype=torch.float32).view(-1, num_q).to(device)
y_valid_tensor = torch.tensor(y_valid_scaled, dtype=torch.float32).view(-1, num_q).to(device)
y_test_tensor  = torch.tensor(y_test_scaled,  dtype=torch.float32).view(-1, num_q).to(device)

# =========================
# Objective of the Bayesian optimisation (variables for each eig_k passed via a closure)
# =========================
def make_objective(num_p_, num_q_, train_loader_, x_valid_tensor_, y_valid_tensor_, device_):
    def objective(trial):
        n_layers      = trial.suggest_int('n_layers', 1, 3)
        hidden_sizes  = [trial.suggest_int(f'n_units_l{i}', 1, 64) for i in range(n_layers)]
        learning_rate = trial.suggest_float('lr', 1e-5, 1e-1, log=True)
        n_epochs      = 100

        set_seed(42)
        model = Net(num_p_, hidden_sizes, num_q_).to(device_)
        optimizer = optim.Adam(model.parameters(), lr=learning_rate)
        criterion = nn.MSELoss()

        model.train()
        for _ in range(n_epochs):
            for batch_x, batch_y in train_loader_:
                optimizer.zero_grad()
                out  = model(batch_x)
                loss = criterion(out, batch_y)
                loss.backward()
                optimizer.step()

        model.eval()
        with torch.no_grad():
            val_loss = criterion(model(x_valid_tensor_), y_valid_tensor_)

        return val_loss.item()
    return objective


optuna.logging.set_verbosity(optuna.logging.WARNING)

# =========================
# Main loop over eig_k
# =========================
summary_rows     = []
scaled_rmse_rows = []   # standardised RMSE (time / SPE / mean) for each eig_k

for eig_k in EIG_K_LIST:
    print("\n" + "=" * 60)
    print(f"Running experiment with eig_k = {eig_k} (top{eig_k} / bottom{eig_k})")
    print("=" * 60)

    foldername = os.path.join(foldername_base, f"eig{eig_k}")
    os.makedirs(foldername, exist_ok=True)

    x_colname = make_x_colname(eig_k)
    print(f"X columns ({len(x_colname)}): {x_colname}")

    x_train = train_df[x_colname]
    x_valid = valid_df[x_colname]
    x_test  = test_df[x_colname]

    num_p = x_train.shape[1]

    # scaler for X fitted separately for each eig_k
    scaler_x = StandardScaler()
    x_train_scaled = scaler_x.fit_transform(x_train)
    x_valid_scaled = scaler_x.transform(x_valid)
    x_test_scaled  = scaler_x.transform(x_test)

    pd.DataFrame({'mean': scaler_x.mean_, 'scale': scaler_x.scale_}).to_csv(
        os.path.join(foldername, 'scaler_x_params_glmnet.csv'), index=False
    )
    pd.DataFrame({'mean': scaler_y.mean_, 'scale': scaler_y.scale_}).to_csv(
        os.path.join(foldername, 'scaler_y_params_glmnet.csv'), index=False
    )

    x_train_tensor = torch.tensor(x_train_scaled, dtype=torch.float32).to(device)
    x_valid_tensor = torch.tensor(x_valid_scaled, dtype=torch.float32).to(device)
    x_test_tensor  = torch.tensor(x_test_scaled,  dtype=torch.float32).to(device)

    train_dataset = TensorDataset(x_train_tensor, y_train_tensor)
    valid_dataset = TensorDataset(x_valid_tensor, y_valid_tensor)
    test_dataset  = TensorDataset(x_test_tensor,  y_test_tensor)

    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)
    valid_loader = DataLoader(valid_dataset, batch_size=batch_size, shuffle=False)
    test_loader  = DataLoader(test_dataset,  batch_size=batch_size, shuffle=False)

    print(f"Train batches: {len(train_loader)}, Valid batches: {len(valid_loader)}, Test batches: {len(test_loader)}")

    # ---- 1) Bayesian optimisation with Optuna (evaluated on VALID) ----
    sampler_EI = optuna.integration.BoTorchSampler(candidates_func=EI_candidates_func)
    db_path    = os.path.join(foldername, f"nas_glmnet_valid_eig{eig_k}.db")

    study = optuna.create_study(
        direction='minimize',
        sampler=sampler_EI,
        storage=f"sqlite:///{db_path}",
        study_name=f"nas_study_valid_eig{eig_k}"
    )

    objective_fn = make_objective(
        num_p, num_q, train_loader,
        x_valid_tensor, y_valid_tensor, device
    )
    study.optimize(objective_fn, n_trials=n_bo_trials)

    trial = study.best_trial
    print(f"Best trial number : {trial.number}")
    print(f"Best Value (VALID MSE): {trial.value:.6f}")
    print("Params:")
    for k, v in trial.params.items():
        print(f"  {k}: {v}")

    # ---- 2) Final training (with early stopping) ----
    best_params   = trial.params
    n_layers      = best_params['n_layers']
    hidden_sizes  = [best_params[f'n_units_l{i}'] for i in range(n_layers)]
    learning_rate = best_params['lr']

    set_seed(seed)
    model_ep  = Net(num_p, hidden_sizes, num_q).to(device)
    optimizer = optim.Adam(model_ep.parameters(), lr=learning_rate)
    criterion = nn.MSELoss()

    best_val_rmse     = float('inf')
    best_epoch        = -1
    no_improve_count  = 0
    best_model_state  = None

    model_ep.train()
    for epoch in range(1, max_epochs + 1):
        for batch_x, batch_y in train_loader:
            optimizer.zero_grad()
            out  = model_ep(batch_x)
            loss = criterion(out, batch_y)
            loss.backward()
            optimizer.step()

        model_ep.eval()
        with torch.no_grad():
            val_pred = model_ep(x_valid_tensor).cpu().numpy()
        val_rmse = np.sqrt(mean_squared_error(y_valid_scaled, val_pred))

        if val_rmse < best_val_rmse:
            best_val_rmse    = val_rmse
            best_epoch       = epoch
            best_model_state = copy.deepcopy(model_ep.state_dict())
            no_improve_count = 0
        else:
            no_improve_count += 1

        model_ep.train()
        if epoch % 50 == 0:
            print(f"Epoch: {epoch:4d}, VALID RMSE (scaled): {val_rmse:.6f}"
                  f"  {'<- best' if epoch == best_epoch else f'(no improve: {no_improve_count}/{patience})'}")

        if no_improve_count >= patience:
            print(f"Early stopping at epoch {epoch} (no improvement for {patience} epochs)")
            break

    print(f"Best epoch (by VALID): {best_epoch}, RMSE: {best_val_rmse:.6f}")

    # ---- 3) Save the best model -> evaluate on TEST ----
    best_model = Net(num_p, hidden_sizes, num_q).to(device)
    best_model.load_state_dict(best_model_state)

    checkpoint = {
        'model_state_dict': best_model.state_dict(),
        'best_params'     : best_params,
        'input_size'      : num_p,
        'output_size'     : num_q,
        'n_layers'        : n_layers,
        'hidden_sizes'    : hidden_sizes,
        'learning_rate'   : learning_rate,
        'n_epochs'        : best_epoch,
        'x_colname'       : x_colname,
        'y_colname'       : y_colname,
        'eig_k'           : eig_k,
    }
    ckpt_path = os.path.join(foldername, f"best_model_glmnet_valid_eig{eig_k}.pth")
    torch.save(checkpoint, ckpt_path)
    print(f"Checkpoint saved: {ckpt_path}")

    best_model.eval()
    with torch.no_grad():
        pred_test = best_model(x_test_tensor).cpu().numpy()

    rmse_test_scaled = np.sqrt(mean_squared_error(y_test_scaled, pred_test))
    print(f"[TEST] RMSE (scaled)  : {rmse_test_scaled:.6f}")

    pred_log  = scaler_y.inverse_transform(pred_test)
    pred_orig = np.exp(pred_log)
    rmse_test_orig = np.sqrt(mean_squared_error(y_test, pred_orig))
    print(f"[TEST] RMSE (original): {rmse_test_orig:.6f}")

    # standardised RMSE of each target from the scaled plot
    y_test_scaled_df = pd.DataFrame(y_test_scaled, columns=y_test.columns)
    rmse_per_target_scaled = plot_predictions_vs_observations(
        y_test_scaled_df, pred_test, folder_name=foldername, name="scaled"
    )
    rmse_per_target_orig = plot_predictions_vs_observations(
        y_test, pred_orig, folder_name=foldername, name="orig"
    )

    # record the standardised RMSE (time / SPE / mean) for this eig_k
    scaled_row = {'eig_k': eig_k}
    for c in y_colname:
        scaled_row[f'{c}_rmse_scaled'] = rmse_per_target_scaled[c]
    scaled_row['mean_rmse_scaled'] = np.mean(
        [rmse_per_target_scaled[c] for c in y_colname]
    )
    scaled_rmse_rows.append(scaled_row)

    summary_rows.append({
        'eig_k'                : eig_k,
        'n_features'           : num_p,
        'best_epoch'           : best_epoch,
        'valid_rmse_scaled'    : best_val_rmse,
        'test_rmse_scaled'     : rmse_test_scaled,
        'test_rmse_orig_time'  : rmse_per_target_orig.get('time_glmnet', np.nan),
        'test_rmse_orig_spe'   : rmse_per_target_orig.get('SPE_glmnet', np.nan),
    })

# =========================
# 4) Save the summary of the eigenvalue-feature comparison
# =========================
summary_df = pd.DataFrame(summary_rows)
summary_path = os.path.join(foldername_base, "eigenvalue_ablation_summary.csv")
summary_df.to_csv(summary_path, index=False)
print("\n" + "=" * 60)
print("Ablation summary:")
print(summary_df)
print(f"Saved summary to: {summary_path}")

# =========================
# 5) Save the standardised RMSE (time / SPE / mean) for each eig_k,
#    sorted by the mean (smallest first)
# =========================
scaled_rmse_df = pd.DataFrame(scaled_rmse_rows).sort_values(
    'mean_rmse_scaled'
).reset_index(drop=True)
scaled_rmse_path = os.path.join(foldername_base, "eigenvalue_scaled_rmse.csv")
scaled_rmse_df.to_csv(scaled_rmse_path, index=False)
print("\n" + "=" * 60)
print("Scaled RMSE (sorted by mean):")
print(scaled_rmse_df)
print(f"Saved scaled RMSE to: {scaled_rmse_path}")