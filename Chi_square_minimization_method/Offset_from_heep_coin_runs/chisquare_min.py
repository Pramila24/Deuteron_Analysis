import numpy as np
import scipy.linalg as la


# Constants
deg_to_rad = np.pi / 180.0
Mp = 0.938272
Eb = 10.6005
Ef = 8.5342488

runs = [3288, 3371, 3374, 3377]
Pp = np.array([2.935545, 3.47588, 2.310387, 1.891227])
th_e = np.array([12.194, 13.930, 9.928, 8.495]) * deg_to_rad
th_p = np.array([37.338, 33.545, 42.9, 47.605]) * deg_to_rad
ph_e = 0
ph_p = 0

# Observables from data and simc
Wdata = np.array([0.942152, 0.940229, 0.943825, 0.955773])
Wdata_err = np.array([9.71583e-05, 3.77469e-04, 6.01623e-05, 4.13486e-05])
Wsimc = np.array([0.94469, 0.94448, 0.943591, 0.942933])
Wsimc_err = np.array([1.38566e-04, 6.68696e-05, 2.81437e-05, 2.23158e-05])

Emdata = np.array([6.58010e-03, 5.28983e-03, 6.26642e-03, 5.32021e-03])
Emdata_err = np.array([4.69566e-05, 1.47900e-04, 3.65940e-05, 4.14000e-05])
Emsimc = np.array([5.51066e-03, 6.06886e-03, 6.58531e-03, 6.47356e-03])
Emsimc_err = np.array([9.05496e-05, 3.84833e-05, 1.74941e-05, 1.37620e-05])

PmXdata = np.array([-1.29097e-03, -9.93560e-04, -5.53961e-04, 7.13300e-03])
PmXdata_err = np.array([3.30069e-05, 1.75464e-04, 2.87979e-05, 2.08787e-05])
PmXsimc = np.array([-5.88228e-04, -5.91870e-04, -6.80114e-04, -1.02935e-03])
PmXsimc_err = np.array([7.19549e-05, 2.48719e-05, 2.31853e-05, 1.11287e-05])

PmZdata = np.array([5.77521e-03, 5.42333e-03, 5.33585e-03, 3.43146e-03])
PmZdata_err = np.array([7.85239e-05, 1.13337e-04, 3.97398e-05, 3.45173e-05])
PmZsimc = np.array([5.82773e-03, 6.50794e-03, 6.35112e-03, 6.02246e-03])
PmZsimc_err = np.array([1.06200e-04, 2.89548e-05, 1.55699e-05, 1.62459e-05])

# Systematic error
p_syst = 1e-4

# Initialize matrices
C = np.zeros((12, 4))
bVec = np.zeros(12)
bVec_err = np.zeros(12)
N_inv = np.zeros((12, 12))

for i in range(3):
    row = 4 * i

    # Observed differences
    dW_obs = Wdata[i] - Wsimc[i]
    dEm_obs = Emdata[i] - Emsimc[i]
    dPmx_obs = PmXdata[i] - PmXsimc[i]
    dPmz_obs = PmZdata[i] - PmZsimc[i]

    # Errors squared
    dW_obs_err2 = Wdata_err[i]**2 + Wsimc_err[i]**2 + (p_syst * Ef)**2
    dEm_obs_err2 = Emdata_err[i]**2 + Emsimc_err[i]**2 + (p_syst * Ef)**2
    dPmx_obs_err2 = PmXdata_err[i]**2 + PmXsimc_err[i]**2 + (p_syst * Ef)**2
    dPmz_obs_err2 = PmZdata_err[i]**2 + PmZsimc_err[i]**2 + (p_syst * Ef)**2

    # Fill coefficient matrix C
    C[row][0] = Ef  # dW/dEb * Eb
    C[row][1] = -Eb  # dW/dEf * Ef
    C[row][2] = -2 * Eb * Ef / Mp * np.sin(th_e[i] / 2) * np.cos(th_e[i] / 2)
    C[row][3] = 0

    C[row+1][0] = Eb
    C[row+1][1] = -Ef
    C[row+1][2] = 0
    C[row+1][3] = 0

    C[row+2][0] = 0
    C[row+2][1] = -np.sin(th_e[i]) * np.cos(th_e[i]) * Ef
    C[row+2][2] = -Ef * np.cos(ph_e) * np.cos(th_e[i])
    C[row+2][3] = -Pp[i] * np.cos(th_p[i]) * np.cos(ph_p)

    C[row+3][0] = Eb
    C[row+3][1] = -np.cos(th_e[i]) * np.cos(ph_e) * Ef
    C[row+3][2] = Ef * np.cos(ph_e) * np.sin(th_e[i])
    C[row+3][3] = Pp[i] * np.sin(th_p[i]) * np.cos(ph_p)

    # Fill bVec and errors
    bVec[row:row+4] = [dW_obs, dEm_obs, dPmx_obs, dPmz_obs]
    bVec_err[row:row+4] = np.sqrt([dW_obs_err2, dEm_obs_err2, dPmx_obs_err2, dPmz_obs_err2])
    N_inv[row:row+4, row:row+4] = np.diag(1.0 / np.array([dW_obs_err2, dEm_obs_err2, dPmx_obs_err2, dPmz_obs_err2]))

C, bVec, bVec_err, N_inv




# SVD solution
U, s, VT = la.svd(C, full_matrices=False)
S_inv = np.diag(1.0 / s)
C_pinv = VT.T @ S_inv @ U.T
aVec = C_pinv @ bVec  # <-- this must run before print()

# Residuals and chi2
residuals = bVec - C @ aVec
chi2 = residuals.T @ N_inv @ residuals
chi2dof = chi2 / (12 - 4)

# Covariance matrix
CT = C.T
R1 = N_inv @ C
V = CT @ R1
Vinv = la.inv(V)
param_errors = np.sqrt(np.diag(Vinv))

# Correlation matrix
Cm = np.zeros_like(Vinv)
for i in range(4):
    for j in range(4):
        Cm[i, j] = Vinv[i, j] / (param_errors[i] * param_errors[j])

# === PRINT RESULTS ===
print("=== Optimized Parameters ===")
print(f"dEb/Eb  = {aVec[0]:.6f} ± {param_errors[0]:.6f}")
print(f"dEf/Ef  = {aVec[1]:.6f} ± {param_errors[1]:.6f}")
print(f"dth_e   = {aVec[2]:.6f} rad ± {param_errors[2]:.6f} rad")
print(f"dth_p   = {aVec[3]:.6f} rad ± {param_errors[3]:.6f} rad")

print("\n=== Chi-squared ===")
print(f"Chi²        = {chi2:.4f}")
print(f"Chi² / dof  = {chi2dof:.4f}")

print("\n=== Correlation Matrix ===")

print(Cm)
print("\n=== Inverse of Observed Covariance Matrix (N_inv) ===")
print(N_inv)

print("\n=== Covariance Matrix of Parameters (Vinv) ===")
print(Vinv)
