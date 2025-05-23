from time import time

import cosmoglobe as cg
import healpy as hp
import matplotlib.pyplot as plt
import numpy as np
from BackgroundCosmology import BackgroundCosmology
from Perturbations import Perturbations
from PowerSpectrum import PowerSpectrum
from RecombinationHistory import RecombinationHistory

t_start = time()

h = 0.67
OmegaB0 = 0.05
OmegaCDM0 = 0.267
OmegaK0 = 0
Neff = 3.046
TCMB = 2.7255

cosmo = BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB)
cosmo.solve()

Yp = 0.245
z_reion = 8
delta_z_reion = 0.5
z_Hereion = 3.5
delta_z_Hereion = 0.5
reion = True

rec = RecombinationHistory(
    cosmo, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, reion
)
rec.solve()

pert = Perturbations(cosmo, rec)
pert.solve()

A_s = 2.1e-9
n_s = 0.965
kpivot_mpc = 0.05

ps = PowerSpectrum(cosmo, rec, pert, A_s, n_s, kpivot_mpc)
ps.solve()

# get eta0
eta0 = cosmo.eta_of_x(0)
H0 = cosmo.get_H_of_x(0)

m = 1.0
s = 1.0
Mpc = 3.08567758e22 * m
k_min = 0.00005 / Mpc
deltak = 2 * np.pi / eta0 / 6
k_max = 3000.0 / eta0
n_k = int((k_max - k_min) / deltak)

c = 2.99792458e8 * m / s

ks = np.linspace(k_min, k_max, n_k)

magic_ells = [6, 100, 200, 500, 1000]
thetaT = np.zeros((len(magic_ells), len(ks)))

plt.figure(figsize=(10, 6))

for elli in range(len(magic_ells)):
    ell = magic_ells[elli]
    for ki in range(len(ks)):
        thetaT[elli, ki] = ps.get_thetaT_ell_of_k(ell, ks[ki])
    plt.plot(ks * c / H0, thetaT[elli, :], label=f"$\ell$ = {ell}")
plt.legend()
plt.xlim(0, 600)
plt.ylim(-0.015, 0.015)
plt.ylabel(r"$\Theta_\ell$")
plt.xlabel(r"$k c / H_0$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/PowerSpectrum/theta_ell_of_k.png"
)
plt.close()

# plot theta^2 / k
plt.figure(figsize=(10, 6))
for elli in range(len(magic_ells)):
    plt.plot(
        c * ks / H0,
        thetaT[elli, :] ** 2 * H0 / (ks * c),
        label=f"$\ell$ = {magic_ells[elli]}",
    )
plt.xlabel(r"$k c / H_0$")
plt.ylabel(r"$\Theta_\ell^2 H_0/ (kc)$")
plt.xlim(0, 600)
plt.ylim(0, 3e-6)
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/PowerSpectrum/theta_ell2_of_k.png"
)
plt.close()

# plot matter power spectrum
plt.figure(figsize=(10, 6))
x = 0
matter_power_spectrum = np.zeros(len(ks))
for ki in range(len(ks)):
    matter_power_spectrum[ki] = ps.get_matter_power_spectrum(x, ks[ki])
matter_power_spectrum = matter_power_spectrum / (Mpc / h) ** 3

# matter-radiation equality
x_eq = -8.1319
a_eq = np.exp(x_eq)
H_eq = cosmo.get_H_of_x(x_eq)
k_eq = a_eq * H_eq / c

plt.plot(ks * Mpc / h, matter_power_spectrum)
plt.vlines(
    k_eq * Mpc / h,
    np.min(matter_power_spectrum),
    np.max(matter_power_spectrum) + 1e4,
    color="red",
    linestyle="--",
    label=r"$k_{eq}$",
)
plt.xlabel(r"$k$ (h/Mpc)")
plt.ylabel(r"$P(k)$ (Mpc/h)$^3$")
plt.xscale("log")
plt.yscale("log")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/PowerSpectrum/matter_power_spectrum.png"
)
plt.close()

# get planck points
file_path = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/planck_lowell_TT.txt"
planck_data = np.loadtxt(file_path, skiprows=1)
planck_ells = planck_data[:, 0]
planck_cell = planck_data[:, 1]
planck_errup = planck_data[:, 2]
planck_errdown = planck_data[:, 3]

file_path = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/planck_highell_TT.txt"
planck_data_high = np.loadtxt(file_path, skiprows=1)
planck_ells = np.append(planck_ells, planck_data_high[:, 0])
planck_cell = np.append(planck_cell, planck_data_high[:, 1])
planck_errup = np.append(planck_errup, planck_data_high[:, 3])
planck_errdown = np.append(planck_errdown, planck_data_high[:, 2])

file_path = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/planck_highell_EE.txt"
planck_data = np.loadtxt(file_path, skiprows=1)
planck_ells_EE = planck_data[:, 0]
planck_cell_EE = planck_data[:, 1]
planck_errdown_EE = planck_data[:, 2]
planck_errup_EE = planck_data[:, 3]

file_path = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/planck_highell_TE.txt"
planck_data = np.loadtxt(file_path, skiprows=1)
planck_ells_TE = planck_data[:, 0]
planck_cell_TE = planck_data[:, 1]
planck_errdown_TE = planck_data[:, 2]
planck_errup_TE = planck_data[:, 3]

ells = np.arange(2, 2000, 1)

# plot power spectrum
plt.figure(figsize=(10, 6))
cell_TT = np.zeros(len(ells))
for elli in range(len(ells)):
    cell_TT[elli] = ps.get_cell_TT(ells[elli])

plt.plot(ells, cell_TT)
plt.errorbar(
    planck_ells,
    planck_cell,
    yerr=[planck_errdown, planck_errup],
    fmt="o",
    markersize=3,
    label="Planck 2018",
)
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\frac{\ell(\ell+1)C_\ell^{TT}}{2\pi}$ ($\mu K^2$)")
plt.xscale("log")
plt.ylim(-1000, 8000)
plt.xlim(1, 3000)
plt.show()

# check factor between curve and points
factor = np.zeros(len(planck_ells))
for elli in range(len(planck_ells)):
    factor[elli] = ps.get_cell_TT(planck_ells[elli]) / planck_cell[elli]
print(f"factor = {np.mean(factor)}")

# plot polarization power spectrum
plt.figure(figsize=(10, 6))
cell_EE = np.zeros(len(ells))
print("ell", "cell_EE")
for elli in range(len(ells)):
    cell_EE[elli] = (
        ps.get_cell_EE(ells[elli])
        #         # * factorial(ells[elli] + 2)
        #         # / factorial(ells[elli] - 2)
    )
#     print(ells[elli], cell_EE[elli])

plt.plot(ells, cell_EE)
plt.errorbar(
    planck_ells_EE,
    planck_cell_EE,
    yerr=[
        planck_errdown_EE,
        planck_errup_EE,
    ],
    fmt="o",
    markersize=3,
    label="Planck 2018",
)
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\frac{\ell(\ell+1)C_\ell^{EE}}{2\pi}$ ($\mu K^2$)")
plt.xlim(0, 2500)
plt.ylim(-20, 120)
plt.legend()
plt.show()

# plot TE power spectrum
cell_TE = np.zeros(len(ells))
for elli in range(len(ells)):
    cell_TE[elli] = ps.get_cell_TE(ells[elli])  # * np.sqrt(
#     # factorial(ells[elli] + 2) / factorial(ells[elli] - 2)
#     # )

plt.figure(figsize=(10, 6))
plt.plot(ells, cell_TE)
plt.errorbar(
    planck_ells_TE,
    planck_cell_TE,
    yerr=[planck_errdown_TE, planck_errup_TE],
    fmt="o",
    markersize=3,
    label="Planck 2018",
)
plt.xlabel(r"$\ell$")
plt.ylabel(r"$\frac{\ell(\ell+1)C_\ell^{TE}}{2\pi}$ ($\mu K^2$)")
plt.xlim(-500, 3000)
plt.ylim(-150, 150)
plt.legend()
plt.show()

# plot maps
alms = hp.sphtfunc.synalm(
    cls=cell_TT / ((ells * (ells + 1)) / (2.0 * np.pi)), lmax=2000
)
maps = hp.sphtfunc.alm2map(alms, nside=512, lmax=2000)
cg.plot(maps, comp="cmb", min=-300, max=200, unit="uK")

alms = hp.sphtfunc.synalm(
    cls=cell_EE / ((ells * (ells + 1)) / (2.0 * np.pi)), lmax=2000
)
maps = hp.sphtfunc.alm2map(alms, nside=512, lmax=2000)
cg.plot(maps, comp="cmb", min=-300, max=200, unit="uK")
plt.show()


t_end = time()
print(f"Total time: {(t_end - t_start)/60} minutes")
