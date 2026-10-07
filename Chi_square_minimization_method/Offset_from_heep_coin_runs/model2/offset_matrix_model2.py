import scipy.linalg as la
#import uproot
import numpy as np
import pandas as pd
from scipy.linalg import pinv
import matplotlib.pyplot as plt
import os
from scipy.stats import norm

#MODEL 1 (EB,thetae) HERE EB IS FIXED AND THETAE IS CHANGING and full derivative of W is taken and the elastic Run momentum is taken
model = "Model 2"

# Constants
deg_to_rad = np.pi / 180.0
Mp = 0.938272
Eb = 10.542  # assumed same for all runs


df = pd.read_csv("/Users/pramila/Research/deut_offline_replay/Optics_study_Chi_square/fit_all_summary_recoil.csv")

# Per-run final energy and angles (from data)
Pf = np.array([3.499,3.145,2.783, 2.418, 2.048,1.664])
th_e = np.array([14.153,12.94, 11.700,10.435, 9.125,7.704]) * deg_to_rad
th_p = np.array([33.344,35.75, 38.549,41.812,45.667,50.500]) * deg_to_rad
# Runs
Run = [20841,20846,20851,20858,20861,20869]
ph_e = 0
ph_p= 0
row = df[df["Runs"] == 20841].iloc[0]
# Systematic error
p_syst = 1e-4


n = len(Run)
W=np.zeros(n)
Ep= np.zeros(n)
# Assuming AX = B + epsilon ; A is the model matrix and X is the vector matrix that needs to be determined and B is the observed quantities
# Want to find the value of epsilon such that epsilon is minimized as such AX-B =0 
# X**2 (Chi-square minimization): (B-AX)^T*N^-1*(B-AX). where N^-1 = 1/N.  ; where N is the sum over (1/sigma squared) is covariance matrix
# Taking W with respect to (Eb,Ef,th_e) where Eb and Ef are constant for all Runs value
A = np.zeros((4*n, 4))
BVector = np.zeros(4*n)
BVector_err = np.zeros(4*n)
N_inv = np.zeros((4*n, 4*n))
Ef= np.array([7.866,8.208,8.55,8.892,9.234,9.576])

for i in range(n):
    row_index = 4*i
    th = th_e[i]
    thp = th_p[i]
    
    # Per row data
    row = df[df["Runs"] == Run[i]].iloc[0]  # select the correct row

    # Calculate Measured_dW and Measured_dW_err2 for this run
    Measured_dW = row["W_data"] - row["W_simc"]
    Measured_dEm = row["Em_data"] - row["Em_simc"]
    Measured_dPmx = row["Pmx_data"] - row["Pmx_simc"]
    Measured_dPmz = row["Pmz_data"] - row["Pmz_simc"]
    Measured_dW_err2 = row["W_data_err"]**2 + row["W_simc_err"]**2 + (p_syst * Ef[i])**2
    Measured_dEm_err2 = row["Em_data_err"]**2 + row["Em_simc_err"]**2 + (p_syst * Ef[i])**2
    Measured_dPmx_err2 = row["Pmx_data_err"]**2 + row["Pmx_simc_err"]**2 + (p_syst * Ef[i])**2
    Measured_dPmz_err2 = row["Pmz_data_err"]**2 + row["Pmz_simc_err"]**2 + (p_syst * Ef[i])**2
    
    arg= Mp*Mp + 2 * Mp* (Eb-Ef[i]) - 4*Eb*Ef[i]*np.sin(th_e[i]/2)**2
    W[i]= np.sqrt(arg)
    Ep[i]= np.sqrt(Mp **2 + Pf[i]**2)
    # Derivatives from W^2 formula
    dW_dEb = ((Mp- 2*Ef[i]*np.sin(th_e[i]/2)**2)/W[i])
    dW_dEf = ((-Mp-2*Eb*np.sin(th_e[i]/2)**2)/W[i])
    dW_dth = ((-4 * Eb * Ef[i] * np.sin(th_e[i]/2) * np.cos(th_e[i]/2)) /W[i])
    dW_dPf = 0
    dEm_dEb = 1
    dEm_dEf= -1
    dEm_dth= 0
    dEm_dPf= -Pf[i]/Ep[i]
    dPmx_dEb= 0
    dPmx_dEf= - np.sin(th)
    dPmx_dth= -Ef[i]*np.cos(th)
    dPmx_dPf= -np.sin(thp)
    dPmz_dEb=1
    dPmz_dEf= -np.cos(th)
    dPmz_dth=Ef[i]*np.sin(th)
    dPmz_dPf=-np.cos(thp)

    #dW_dEb = (Ef / Eb)
    # dW_dEf = (-Eb / Ef[i])
    # dW_dth = (-2 * Eb * Ef[i] * np.sin(th)) / (Mp)

    A[row_index][0] = dW_dEb * Eb
    A[row_index][1] = dW_dEf * Ef[i]
    A[row_index][2] = dW_dth
    A[row_index][3] = dW_dPf *Pf[i]
    A[row_index+1][0] = dEm_dEb * Eb
    A[row_index+1][1] = dEm_dEf * Ef[i]
    A[row_index+1][2] = dEm_dth
    A[row_index+1][3] = dEm_dPf *Pf[i]
    A[row_index+2][0] = dPmx_dEb * Eb
    A[row_index+2][1] = dPmx_dEf * Ef[i]
    A[row_index+2][2] = dPmx_dth
    A[row_index+2][3] = dPmx_dPf *Pf[i]
    A[row_index+3][0] = dPmz_dEb * Eb
    A[row_index+3][1] = dPmz_dEf * Ef[i]
    A[row_index+3][2] = dPmz_dth
    A[row_index+3][3] = dPmz_dPf *Pf[i]

    BVector[row_index:row_index+4] = [Measured_dW,Measured_dEm,Measured_dPmx,Measured_dPmz]
    BVector_err[row_index:row_index+4] = np.sqrt([Measured_dW_err2,Measured_dEm_err2,Measured_dPmx_err2,Measured_dPmz_err2] )
    N_inv[row_index:row_index+4 , row_index:row_index+4] = np.diag(1.0 / np.array([Measured_dW_err2,Measured_dEm_err2,Measured_dPmx_err2,Measured_dPmz_err2]))


# U, s, VT = la.svd(A, full_matrices=False)
U, s, VT = la.svd(A, full_matrices=False)


N_sqrt = np.sqrt(N_inv)  # Element-wise sqrt for diagonal N_inv

# But N_inv is a 4n x 4n matrix, so W_sqrt is also 4n x 4n diagonal matrix
# So we do matrix multiplication
A_weighted = N_sqrt @ A
B_weighted = N_sqrt @ BVector

# # Then SVD of weighted A
U_w, s_w, VT_w = la.svd(A_weighted, full_matrices=False)
S_inv_w = np.diag(1.0 / s_w)
C_pinv_w = VT_w.T @ S_inv_w @ U_w.T

XVector_svd_weighted = C_pinv_w @ B_weighted


S_inv = np.diag(1.0 / s)
C_pinv = VT.T @ S_inv @ U.T
XVector_svd = C_pinv @ BVector
XVector_normal = np.linalg.inv(A.T @ N_inv @ A) @ (A.T @ N_inv) @ BVector

# Residuals and chi²
# Assuming AX = B + epsilon ; A is the model matrix and X is the vector matrix that needs to be determined and B is the observed quantities
# Want to find the value of epsilon such that epsilon is minimized as such AX-B =0 
# X**2 (Chi-square minimization): (B-AX)^T*N^-1*(B-AX). where N^-1 = 1/N.  ; where N is the sum over (1/sigma squared) is covariance matrix
residuals_svd = BVector - A @ XVector_svd
chi2_svd = residuals_svd.T @ N_inv @ residuals_svd
chi2dof_svd = chi2_svd / (4*n - 4)

residuals_normal = BVector - A @ XVector_normal
chi2_normal = residuals_svd.T @ N_inv @ residuals_normal
chi2dof_normal = chi2_normal / (4*n - 4)

residuals_svd_weighted = BVector - A @ XVector_svd_weighted
chi2_svd_weighted = residuals_svd.T @ N_inv @ residuals_svd_weighted
chi2dof_svd_weighted = chi2_svd_weighted / (4*n - 4)

# Covariance matrix.      (Transpose of A is A.T ). We have H=(A^T*W*A) where W= N^-1
AT = A.T
H = AT @ N_inv @ A
Hinv = np.linalg.pinv(H)
V= Hinv
param_errors = np.sqrt(np.diag(Hinv))  # V[a]= Hinv is the Covariance matrix a by "error" propagation, here H is Hessian matrix defined as H= d^2F/dai.daj.  F= 1/2*S(a) S(a)= (B-AX)^T*W*(B-AX)

# Correlation matrix
Cm = np.zeros_like(Hinv)
for i in range(4):
    for j in range(4):
        Cm[i, j] = Hinv[i, j] / (param_errors[i] * param_errors[j])
print()
print(A)
print()
print("Matrix of Measured Parameters:BVector")
print(BVector)
print()
print("Condition number of A:", np.linalg.cond(A))
print()
print("Singular values:", s)
print()
print("====W====")
print(f"W= {W}")

# === Print Results ===
print()

#Print both
# print("XVector using SVD:")
# print(XVector_svd_weighted)

# print("XVector using SVD without weight:")
# print(XVector_svd)

# print("\nXVector using Normal Equations:")
# print(XVector_normal)

# # Optional: print the difference between the two
# print("\nDifference (SVD - Normal):")
# print(XVector_svd - XVector_normal)
print("=== Optimized Parameters Using SVD without weight===")
print(f"dEb/Eb  = {XVector_svd[0]:.6f} ± {param_errors[0]:.6f}")
print(f"dEf/Ef  = {XVector_svd[1]:.6f} ± {param_errors[1]:.6f}")
print(f"dθ_e    = {XVector_svd[2]:.6f} rad ± {param_errors[2]:.6f} rad")
print(f"dPf/Pf  = {XVector_svd[3]:.6f} ± {param_errors[3]:.6f}")
print()
print("=== Optimized Parameters Using the definition  ===")
print(f"dEb/Eb  = {XVector_normal[0]:.6f} ± {param_errors[0]:.6f}")
print(f"dEf/Ef  = {XVector_normal[1]:.6f} ± {param_errors[1]:.6f}")
print(f"dθ_e    = {XVector_normal[2]:.6f} rad ± {param_errors[2]:.6f} rad")
print(f"dPf/Pf  = {XVector_normal[3]:.6f} ± {param_errors[3]:.6f}")
print()
print("=== Optimized Parameters Using SVD  applying weight===")
print(f"dEb/Eb  = {XVector_svd_weighted[0]:.6f} ± {param_errors[0]:.6f}")
print(f"dEf/Ef  = {XVector_svd_weighted[1]:.6f} ± {param_errors[1]:.6f}")
print(f"dθ_e    = {XVector_svd_weighted[2]:.6f} rad ± {param_errors[2]:.6f} rad")
print(f"dPf/Pf  = {XVector_svd_weighted[3]:.6f} ± {param_errors[3]:.6f}")
print()

print("\n=== Chi-squared normal ===")
print(f"Chi²        = {chi2_normal:.4f}")
print(f"Chi² / dof  = {chi2dof_normal:.4f}")
print()
print("\n=== Chi-squared svd ===")
print(f"Chi²        = {chi2_svd:.4f}")
print(f"Chi² / dof  = {chi2dof_svd:.4f}")
print()
print("\n=== Chi-squared svd With weight ===")
print(f"Chi²        = {chi2_svd_weighted:.4f}")
print(f"Chi² / dof  = {chi2dof_svd_weighted:.4f}")
print()

plt.plot(BVector -A @ XVector_svd, "r--" ,BVector -A @ XVector_svd_weighted , "b", BVector- A @ XVector_normal, "g")
plt.legend()
plt.show()
print("\n=== Correlation Matrix ===")
print(Cm)
print()
print("\n=== Inverse of Observed Covariance Matrix (N_inv) ===")
print(N_inv)
print()
print("\n=== Covariance Matrix of Parameters (Vinv) ===")
print(Hinv)
