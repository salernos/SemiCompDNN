#===============================================================================
#
#  PROGRAM: METHOD_DNN.py
#
#  AUTHORS: Stephen Salerno, Zhilin Zhang
#
#  PURPOSE: Generate simulation results for the proposed method
#
#  INPUT:   simdat{X}.RData
#
#  OUTPUT:  Sim_05_{i}_{sim}.pkl
#
#  NOTES:
#
#  The simulation results are stored in 16 x nsims .pkl files, with each file 
#  corresponding to one replicate from one simulation setting. Each file is an 
#  n x 10 array, storing the estimated frailty variance, baseline hazard 
#  parameters, and log risk functions for each simulated dataset.
#
#===============================================================================

#=== INITIALIZATION ============================================================

#--- NECESSARY PACKAGES --------------------------------------------------------

import torch
import torch.nn as nn
import torch.optim as optim
from torch.utils.data import Dataset, DataLoader
import numpy as np
import pandas as pd
import math
import matplotlib.pyplot as plt
from scipy.special import psi, gammaln
from scipy.special import digamma
import pyreadr
import os
import time
from tqdm import tqdm
import pickle

#--- SETUP SLURM ---------------------------------------------------------------

index = int(os.getenv('SLURM_ARRAY_TASK_ID'))

#=== HELPER FUNCTIONS ==========================================================

#--- LOSS FUNCTIONS ------------------------------------------------------------

def loss1_fun_stable(h1, d1, y1, gamma):
    h1 = h1.view(-1)
    d1 = d1.view(-1)
    y1 = y1.view(-1)
    
    gamma = gamma.view(-1)

    risk_matrix = (y1.view(-1, 1) <= y1.view(1, -1)).float()

    h1_stable = h1 - torch.max(h1)
    
    exp_h1 = gamma * torch.exp(h1_stable)

    denom = torch.matmul(risk_matrix, exp_h1) + 1e-8
    
    log_risk = torch.log(denom) + torch.max(h1)

    loss_vector = -d1 * (h1 - log_risk)

    return loss_vector.sum() / (d1.sum() + 1e-8)

def loss2_fun_stable(h2, d1, d2, y1, gamma):
    h2 = h2.view(-1)
    d1 = d1.view(-1)
    d2 = d2.view(-1)
    y1 = y1.view(-1)
    
    gamma = gamma.view(-1)

    event_mask = ((d1 == 0) & (d2 == 1)).float()
    
    risk_matrix = (y1.view(-1, 1) <= y1.view(1, -1)).float()

    h2_stable = h2 - torch.max(h2)
    
    exp_h2 = gamma * torch.exp(h2_stable)

    denom = torch.matmul(risk_matrix, exp_h2) + 1e-8
    
    log_risk = torch.log(denom) + torch.max(h2)

    loss_vector = -event_mask * (h2 - log_risk)

    return loss_vector.sum() / (event_mask.sum() + 1e-8)

def loss3_fun_stable(h3, d1, d2, y1, y2, gamma):
    h3 = h3.view(-1)
    d1 = d1.view(-1)
    d2 = d2.view(-1)
    y1 = y1.view(-1)
    y2 = y2.view(-1)
    
    gamma = gamma.view(-1)

    event_mask = ((d1 == 1) & (d2 == 1)).float()

    risk_matrix = ((y1.view(1, -1) < y2.view(-1, 1)) &
                   (y2.view(1, -1) >= y2.view(-1, 1))).float()

    h3_stable = h3 - torch.max(h3)
    
    exp_h3 = gamma * d1 * torch.exp(h3_stable)

    denom = torch.matmul(risk_matrix, exp_h3) + 1e-8
    
    log_risk = torch.log(denom) + torch.max(h3)

    loss_vector = -event_mask * (h3 - log_risk)

    return loss_vector.sum() / (event_mask.sum() + 1e-8)

#--- MAIN FUNCTION -------------------------------------------------------------

def fit_dnn(formula, 
            data, 
            na_action="na.fail", 
            subset=None,
            dim_layers=[128, 64, 16], 
            lr=0.01, 
            dr=0.1,
            max_epochs=250, 
            max_epochs_theta=100, 
            max_epochs_n=5, 
            verbose=True,
            ll=1, tol=1e-6, 
            theta0=0.5, 
            lr_theta=0.01, 
            batch_size=128):
    #--- Assertions ------------------------------------------------------------
    
    if na_action not in ["na.fail", "na.omit"]:
        raise ValueError('na_action should be either "na.fail" or "na.omit"')
      
    #--- SETUP -----------------------------------------------------------------

    #-- Pre-Process Data
    
    # Assume data are in a pandas DataFrame with the following columns

    true_h1 = data['h1'].values
    true_h2 = data['h2'].values
    true_h3 = data['h3'].values
    
    # Outcomes
    
    y1 = data['Y1'].values
    d1 = data['D1'].values
    y2 = data['Y2'].values
    d2 = data['D2'].values

    # Unique Failure Times by Transition
    
    t1_obs = y1[d1 == 1]
    t2_obs = y2[(d1 == 0) & (d2 == 1)]
    t3_obs = y2[(d1 == 1) & (d2 == 1)]
    
    t1 = np.unique(np.sort(t1_obs))
    t2 = np.unique(np.sort(t2_obs))
    t3 = np.unique(np.sort(t3_obs))
    
    tol_val = tol

    # Features: Same Predictors for All Three Transitions
    
    X1_mat = data[form['X']].values
    X2_mat = data[form['X']].values
    X3_mat = data[form['X']].values
    
    n = X1_mat.shape[0]

    #-- Initialize Parameters
    
    theta = torch.tensor(theta0, requires_grad=True, dtype=torch.float32)
    
    optimizer_theta = optim.Adam([theta], lr=lr_theta)

    # Baseline Hazards (Nelson-Aalen)
    
    d1_j = np.array([np.sum(d1 * (y1 == t)) for t in t1])
    n1_j = np.array([np.sum(y1 >= t) for t in t1])
    
    lam01 = d1_j / n1_j

    d2_j = np.array([np.sum((1 - d1) * d2 * (y2 == t)) for t in t2])
    n2_j = np.array([np.sum((y2 >= t) & (y1 >= t)) for t in t2])
    
    lam02 = d2_j / n2_j

    d3_j = np.array([np.sum(d1 * d2 * (y2 == t)) for t in t3])
    n3_j = np.array([np.sum(y2[d1==1] >= t) for t in t3])
    
    lam03 = d3_j / n3_j

    #-- Neural Network Sub-Architectures
    
    input_dim = X1_mat.shape[1]
    
    def build_model(input_dim, dim_layers, dropout):
        return nn.Sequential(
            nn.Linear(input_dim, dim_layers[0], bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(dim_layers[0], dim_layers[1], bias=False),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(dim_layers[1], dim_layers[2], bias=False),
            nn.Linear(dim_layers[2], 1, bias=False)
        )
    
    h1_model = build_model(input_dim, dim_layers, dr)
    h2_model = build_model(input_dim, dim_layers, dr)
    h3_model = build_model(input_dim, dim_layers, dr)
    
    optimizer_h1 = optim.Adam(h1_model.parameters(), lr=lr)
    optimizer_h2 = optim.Adam(h2_model.parameters(), lr=lr)
    optimizer_h3 = optim.Adam(h3_model.parameters(), lr=lr)

    # Initial Risk Function Values
    
    h1 = h1_model(torch.tensor(X1_mat, dtype=torch.float32))
    h2 = h2_model(torch.tensor(X2_mat, dtype=torch.float32))
    h3 = h3_model(torch.tensor(X3_mat, dtype=torch.float32))

    #-- Internal Helper Functions for Baseline Hazards
    
    def lam01_fun(t, gamma_val, h1_vals):
        numer = np.sum(d1 * (y1 == t))
        h1_np = h1_vals.detach().cpu().numpy().flatten()
        denom = (np.sum(gamma_val * 
          ((y1 >= t).astype(float)) * np.exp(h1_np)))
        return numer/denom if numer > 0 else 0

    def lam02_fun(t, gamma_val, h2_vals):
        numer = np.sum((1 - d1) * d2 * (y1 == t))
        h2_np = h2_vals.detach().cpu().numpy().flatten()
        denom = (np.sum(gamma_val * 
          (((y2 >= t) & (y1 >= t)).astype(float)) * np.exp(h2_np)))
        return numer/denom if numer > 0 else 0

    def lam03_fun(t, gamma_val, h3_vals):
        numer = np.sum(d1 * d2 * (y2 == t))
        h3_np = h3_vals.detach().cpu().numpy().flatten()
        denom = (np.sum(gamma_val * d1 * 
          ((y2 >= t).astype(float) - (y1 >= t).astype(float)) * np.exp(h3_np)))
        return numer/denom if numer > 0 else 0
    
    def Lam01_fun(t, lam01):
        return np.sum(lam01[t1 - t < tol_val])
    def Lam02_fun(t, lam02):
        return np.sum(lam02[t2 - t < tol_val])
    def Lam03_fun(t, lam03):
        return np.sum(lam03[t3 - t < tol_val])

    #-- Theta Contribution to Expected Log-Likelihood
    
    def loss4(theta_tensor, gamma_tensor, log_gamma_tensor):
        loss = (-1 / theta_tensor * 
          torch.log(theta_tensor) + (1/theta_tensor - 1) * 
          log_gamma_tensor - 1/theta_tensor * 
          gamma_tensor - torch.lgamma(1/theta_tensor))
        return -torch.sum(loss)

    #-- Dataset Class for Dataloader
    
    class CustomDataset(Dataset):
        def __init__(self, inputs, outputs):
            self.inputs = inputs
            self.outputs = outputs
        def __len__(self):
            return self.inputs.shape[0]
        def __getitem__(self, idx):
            return self.inputs[idx], self.outputs[idx]

    #--- Neural EM-Algorithm Loop ----------------------------------------------

    # Initialize Loss Tracking
    
    loss_h1_epoch = [np.inf]*(max_epochs+1)
    loss_h2_epoch = [np.inf]*(max_epochs+1)
    loss_h3_epoch = [np.inf]*(max_epochs+1)
    
    loss_theta_epoch = [np.inf]*(max_epochs+1)
    
    loss_theta_inc = 0
    
    loss_h1_inc = 0
    loss_h2_inc = 0
    loss_h3_inc = 0
    
    # Gamma Starting Values
    
    gamma_val = np.ones_like(y1)
    
    log_gamma_val = np.zeros_like(y1)
    
    diff_EM = 100
    
    epoch = 2

    while epoch < max_epochs and diff_EM > tol_val:
      
        #-- E-STEP
        
        theta_val = theta.item()
        
        # Update Risk Functions
        
        h1_all = h1_model(torch.tensor(X1_mat, dtype=torch.float32)).detach()
        h2_all = h2_model(torch.tensor(X2_mat, dtype=torch.float32)).detach()
        h3_all = h3_model(torch.tensor(X3_mat, dtype=torch.float32)).detach()
        
        # Compute Cumulative Hazards
        
        Lam01 = np.array([Lam01_fun(t, lam01) for t in y1])
        Lam02 = np.array([Lam02_fun(t, lam02) for t in y1])
        
        Lam03_y2 = np.array([Lam03_fun(t, lam03) for t in y2])
        Lam03_y1 = np.array([Lam03_fun(t, lam03) for t in y1])
        
        h1_np = h1_all.numpy().flatten()
        h2_np = h2_all.numpy().flatten()
        h3_np = h3_all.numpy().flatten()
        
        # Update Gamma Parameters
        
        a = 1/theta_val + d1 + d2 
        
        b = (1/theta_val + Lam01 * np.exp(h1_np) + 
          Lam02 * np.exp(h2_np) + 
          d1 * (Lam03_y2 - Lam03_y1) * np.exp(h3_np))
          
        gamma_val = a / b
        
        log_gamma_val = digamma(a) - np.log(b)

        print("mean:", np.mean(gamma_val))
        print("var:", np.var(gamma_val))
        
        #-- M-STEP

        # Update Baseline Hazards
        
        lam01_new = np.array([lam01_fun(t, gamma_val, h1_all) for t in t1])
        lam02_new = np.array([lam02_fun(t, gamma_val, h2_all) for t in t2])
        lam03_new = np.array([lam03_fun(t, gamma_val, h3_all) for t in t3])
        
        diff_EM = max(np.max(np.abs((lam01 - lam01_new) / (lam01 + 1e-8))),
                      np.max(np.abs((lam02 - lam02_new) / (lam02 + 1e-8))),
                      np.max(np.abs((lam03 - lam03_new) / (lam03 + 1e-8))))
                      
        lam01 = lam01_new
        lam02 = lam02_new
        lam03 = lam03_new
        
        #-- N-STEP

        # Update Frailty Variance
        
        loss_theta_inc_inner = 0
        
        prev_loss_theta_inner = float('inf')
        
        if loss_theta_inc < 2:
            for epoch_theta in range(max_epochs_theta):
                optimizer_theta.zero_grad()
                loss_theta = loss4(theta,
                  torch.tensor(gamma_val, dtype=torch.float32),
                  torch.tensor(log_gamma_val, dtype=torch.float32))
                loss_theta.backward(retain_graph=True)
                optimizer_theta.step()
                current_loss_theta = loss_theta.item()
                if abs(current_loss_theta - prev_loss_theta_inner) <= 1e-4:
                    break
                
                elif current_loss_theta > prev_loss_theta_inner:
                    loss_theta_inc_inner += 1
                    if loss_theta_inc_inner > 2:
                        break
                if verbose:
                    print(f"Theta Epoch {epoch_theta}")
                    print(f"Loss: {current_loss_theta}")
                    print(f"Theta: {theta.item()}")
                
                prev_loss_theta_inner = current_loss_theta
            
            if current_loss_theta > loss_theta_epoch[epoch-1]:
                loss_theta_inc += 1
            else:
                loss_theta_inc = 0
            
            loss_theta_epoch[epoch] = current_loss_theta

        # Update Neural Network Parameters
        
        # Datasets for Each Transition with Outcomes, Gamma, True h
        
        outputs1 = np.column_stack((y1, d1, y2, d2, gamma_val, true_h1))
        dataset1 = CustomDataset(torch.tensor(X1_mat, dtype=torch.float32),
                                 torch.tensor(outputs1, dtype=torch.float32))
        dl1 = DataLoader(dataset1, 
          batch_size=batch_size, shuffle=True, drop_last=True)

        outputs2 = np.column_stack((y1, d1, y2, d2, gamma_val, true_h2))
        dataset2 = CustomDataset(torch.tensor(X2_mat, dtype=torch.float32),
                                 torch.tensor(outputs2, dtype=torch.float32))
        dl2 = DataLoader(dataset2, 
          batch_size=batch_size, shuffle=True, drop_last=True)

        outputs3 = np.column_stack((y1, d1, y2, d2, gamma_val, true_h3))
        dataset3 = CustomDataset(torch.tensor(X3_mat, dtype=torch.float32),
                                 torch.tensor(outputs3, dtype=torch.float32))
        dl3 = DataLoader(dataset3, 
          batch_size=batch_size, shuffle=True, drop_last=True)
        
        loss1_val = float('inf')
        loss2_val = float('inf')
        loss3_val = float('inf')
        
        if loss_h1_inc < 2:
            
            loss_h1_epoch_n = 0

            num_batches1 = 0

            # Transition 1
            
            for batch_X1, batch_out in dl1:
                num_batches1 += 1
                batch_y1 = batch_out[:, 0]
                batch_d1 = batch_out[:, 1]
                batch_y2 = batch_out[:, 2]
                batch_d2 = batch_out[:, 3]
                batch_gamma = batch_out[:, 4]
                batch_true_h1 = batch_out[:, 5]
                h1_batch = h1_model(batch_X1)
                event_idx = (batch_d1 == 1).nonzero(as_tuple=True)[0]
                loss_h1 = 0.0
                if len(event_idx) > 0:
                    loss_h1 = loss1_fun_stable(
                      h1_batch, batch_d1, batch_y1, batch_gamma)
                    optimizer_h1.zero_grad()
                    loss_h1.backward(retain_graph=True)
                    optimizer_h1.step()

                loss_h1_epoch_n += loss_h1.item()

            if num_batches1 > 0:
                loss1_val = loss_h1_epoch_n / num_batches1
            
            if loss1_val > loss_h1_epoch[epoch-1]:
                loss_h1_inc += 1
            else:
                loss_h1_inc = 0

            loss_h1_epoch[epoch] = loss1_val
            
            h1 = h1_model(torch.tensor(X1_mat, dtype=torch.float32))
        
        if loss_h2_inc < 2:
            
            loss_h2_epoch_n = 0
            
            num_batches2 = 0
            
            # Transition 2
            
            for batch_X2, batch_out in dl2:
                num_batches2 += 1
                batch_y1 = batch_out[:, 0]
                batch_d1 = batch_out[:, 1]
                batch_y2 = batch_out[:, 2]
                batch_d2 = batch_out[:, 3]
                batch_gamma = batch_out[:, 4]
                batch_true_h2 = batch_out[:, 5]
                h2_batch = h2_model(batch_X2)
                event_idx = ((batch_d1 == 0) & 
                             (batch_d2 == 1)).nonzero(as_tuple=True)[0]
                loss_h2 = 0.0
                if len(event_idx) > 0:
                    loss_h2 = loss2_fun_stable(
                      h2_batch, batch_d1, batch_d2, batch_y1, batch_gamma)
                    optimizer_h2.zero_grad()
                    loss_h2.backward(retain_graph=True)
                    optimizer_h2.step()
                loss_h2_epoch_n += loss_h2.item()
                
            if num_batches2 > 0:
                loss2_val = loss_h2_epoch_n / num_batches2
                
            if loss2_val > loss_h2_epoch[epoch-1]:
                loss_h2_inc += 1
            else:
                loss_h2_inc = 0
            
            loss_h2_epoch[epoch] = loss2_val
            
            h2 = h2_model(torch.tensor(X2_mat, dtype=torch.float32))
            
        if loss_h3_inc < 2:
            
            loss_h3_epoch_n = 0
            
            num_batches3 = 0
            
            # Transition 3
            
            for batch_X3, batch_out in dl3:
                num_batches3 += 1
                batch_y1 = batch_out[:, 0]
                batch_d1 = batch_out[:, 1]
                batch_y2 = batch_out[:, 2]
                batch_d2 = batch_out[:, 3]
                batch_gamma = batch_out[:, 4]
                batch_true_h3 = batch_out[:, 5]
                h3_batch = h3_model(batch_X3)
                event_idx = ((batch_d1 == 1) & 
                             (batch_d2 == 1)).nonzero(as_tuple=True)[0]
                loss_h3 = 0.0
                if len(event_idx) > 0:
                    loss_h3 = loss3_fun_stable(
                      h3_batch, batch_d1, batch_d2, 
                      batch_y1, batch_y2, batch_gamma)
                    optimizer_h3.zero_grad()
                    loss_h3.backward(retain_graph=True)
                    optimizer_h3.step()
                loss_h3_epoch_n += loss_h3.item()
            
            if num_batches3 > 0:
                loss3_val = loss_h3_epoch_n / num_batches3
                
            if loss3_val > loss_h3_epoch[epoch-1]:
                loss_h3_inc += 1
            else:
                loss_h3_inc = 0

            loss_h3_epoch[epoch] = loss3_val
            h3 = h3_model(torch.tensor(X3_mat, dtype=torch.float32))

        # Plot Baseline Hazards
        
        if verbose:
          plt.figure(figsize=(15, 4))
          plt.subplot(1, 3, 1)
          plt.plot(t1, np.cumsum(lam01), marker='o')
          plt.title('Baseline Hazards 1')
          plt.xlabel('Failure Time')
          plt.ylabel('Cumulative Baseline Hazards')
          plt.axline((0, 0), slope=2, linestyle='--')
          plt.subplot(1, 3, 2)
          plt.plot(t2, np.cumsum(lam02), marker='o')
          plt.title('Baseline Hazards 2')
          plt.xlabel('Failure Time')
          plt.ylabel('Cumulative Baseline Hazards')
          plt.axline((0, 0), slope=3, linestyle='--')
          plt.subplot(1, 3, 3)
          plt.plot(t3, np.cumsum(lam03), marker='o')
          plt.title('Baseline Hazards 3')
          plt.xlabel('Failure Time')
          plt.ylabel('Cumulative Baseline Hazards')
          plt.axline((0, 0), slope=2, linestyle='--')
          plt.show()
  
          print(f"Epoch {epoch}")
          print(f"loss1: {loss1_val}")
          print(f"loss2: {loss2_val}")
          print(f"loss3: {loss3_val}")
        
        epoch += 1
        
        diff_EM = 100  # Reset diff_EM for next outer iteration

        epochs = range(1, len(loss_h1_epoch))
        
        plt.figure(figsize=(12, 4))

        plt.plot(epochs, loss_h1_epoch[1:], label="Loss h1")
        plt.plot(epochs, loss_h2_epoch[1:], label="Loss h2")
        plt.plot(epochs, loss_h3_epoch[1:], label="Loss h3")
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.title("Loss Curves")
        plt.legend()
        plt.tight_layout()
        plt.show()
        
        if loss_h1_inc >= 2 and loss_h2_inc >= 2 and loss_h3_inc >=2:
            
            break

    return {
        "theta": theta.item(),
        "lam01": lam01,
        "lam02": lam02,
        "lam03": lam03,
        "h1": h1,
        "h2": h2,
        "h3": h3,
        "gamma": gamma_val
    }
    
#=== SIMULATION FUNCTION =======================================================

def run_simulation(sim, dat, form):
    dat_sim = pd.DataFrame(dat[:, :, sim])
    dat_sim.columns = ["Y1", "D1", "Y2", "D2"] + \
                      [f"X{i}" for i in range(1, 13)] + \
                      ["h1", "h2", "h3", "hc", "gamma"] + \
                      ["Ti1", "Ti2", "Ti3", "ctime", "cind"] + \ 
                      ["setting", "n", "theta", "risk", "cens"]
    
    fit_00 = fit_dnn(form, dat_sim,
                     na_action="na.fail", subset=None, 
                     dim_layers=[256, 128, 16],
                     lr=0.0001, dr=0.3,
                     max_epochs=5, max_epochs_theta=10, max_epochs_n=10,
                     batch_size=500, verbose=True, ll=1, tol=1e-6,
                     theta0=1, lr_theta=0.05)
    
    result = {
    "theta": fit_00["theta"],
    "lam01": fit_00["lam01"],
    "lam02": fit_00["lam02"],
    "lam03": fit_00["lam03"],
    "h1": fit_00["h1"].detach().numpy(),
    "h2": fit_00["h2"].detach().numpy(),
    "h3": fit_00["h3"].detach().numpy(),
    "gamma": fit_00["gamma"]
}
    return result

if __name__ == '__main__':
    
    i = 1
    
    print(f"\nSetting {i} ---\n")

    path = os.path.join("Data",f"simdat{i}.RData")
    
    dat = pyreadr.read_r(path)
    dat = dat['dat'].values
    
    form = {
    'Y': ["Y1", "D1", "Y2", "D2"],
    'X': [f"X{i}" for i in range(1, 13)]
    }
    
    nsims = 500  # Number of Simulations

    results = []

    start_time = time.time()
    
    sim = index

    res = run_simulation(sim, dat, form)
    
    os.makedirs("results", exist_ok=True)
    
    with open(os.path.join("results", f"Sim_05_{i}_{sim}.pkl"), "wb") as f:
        pickle.dump(res, f)

    print("Time taken:", time.time() - start_time)
    
#=== END =======================================================================
