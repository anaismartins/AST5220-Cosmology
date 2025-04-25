import matplotlib.pyplot as plt
import numpy as np
from BackgroundCosmology import BackgroundCosmology
from matplotlib import rc
from Perturbations import Perturbations
from RecombinationHistory import RecombinationHistory

m = 1.0
Mpc = 3.08567758e22 * m

k = [0.001 / Mpc, 0.01 / Mpc, 0.1 / Mpc]

h = 0.67
OmegaB0 = 0.05
OmegaCDM0 = 0.267
OmegaK0 = 0
Neff = 3.046
TCMB = 2.7255

Yp = 0.245
z_reion = 8
delta_z_reion = 0.5
z_Hereion = 3.5
delta_z_Hereion = 0.5

x_start = -20
x_end = 0
npts = 5000

cosmo = BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB)
cosmo.solve()

rec = RecombinationHistory(
    cosmo=cosmo,
    Yp=Yp,
    z_reion=z_reion,
    delta_z_reion=delta_z_reion,
    z_Hereion=z_Hereion,
    delta_z_Hereion=delta_z_Hereion,
    reion=True,
)
rec.solve()

pert = Perturbations(
    cosmo=cosmo,
    rec=rec,
)
pert.solve()

delta_gamma = np.zeros((3, npts))
delta_cdm = np.zeros((3, npts))
delta_b = np.zeros((3, npts))
delta_nu = np.zeros((3, npts))

x = np.linspace(x_start, x_end, npts)

for ki in range(len(k)):
    for xi in range(npts):
        delta_gamma[ki, xi] = 4 * pert.get_Theta(x[xi], k[ki], 0)
        delta_cdm[ki, xi] = pert.get_delta_cdm(x[xi], k[ki])
        delta_b[ki, xi] = pert.get_delta_b(x[xi], k[ki])
        delta_nu[ki, xi] = 4 * pert.get_Nu(x[xi], k[ki], 0)

rc("text", usetex=True)
rc("font", family="serif")

# plot delta_gamma
plt.figure(figsize=(8, 6))
plt.plot(x, delta_gamma[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, delta_gamma[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, delta_gamma[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\delta_\gamma$")
plt.title(r"$\delta_\gamma$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/delta_gamma.png"
)

# plot delta_cdm and delta_b together
plt.figure(figsize=(8, 6))
plt.plot(x, delta_cdm[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$", color="C0")
plt.plot(x, delta_b[0], color="C0", linestyle="--")
plt.plot(x, delta_cdm[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$", color="C1")
plt.plot(x, delta_b[1], color="C1", linestyle="--")
plt.plot(x, delta_cdm[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$", color="C2")
plt.plot(x, delta_b[2], color="C2", linestyle="--")
plt.xlabel(r"$x$")
plt.ylabel(r"$\delta$")
plt.yscale("log")
plt.title(r"$\delta_{CDM}$ and $\delta_{b}$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/delta_cdm_delta_b.png"
)

# plot delta_nu
plt.figure(figsize=(8, 6))
plt.plot(x, delta_nu[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, delta_nu[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, delta_nu[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\delta_\nu$")
plt.title(r"$\delta_\nu$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/delta_nu.png"
)

# get velocity perturbations
v_gamma = np.zeros((3, npts))
v_cdm = np.zeros((3, npts))
v_b = np.zeros((3, npts))
v_nu = np.zeros((3, npts))

for ki in range(len(k)):
    for xi in range(npts):
        v_gamma[ki, xi] = -3 * pert.get_Theta(x[xi], k[ki], 1)
        v_cdm[ki, xi] = pert.get_v_cdm(x[xi], k[ki])
        v_b[ki, xi] = pert.get_v_b(x[xi], k[ki])
        v_nu[ki, xi] = -3 * pert.get_Nu(x[xi], k[ki], 1)

# plot v_gamma
plt.figure(figsize=(8, 6))
plt.plot(x, v_gamma[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, v_gamma[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, v_gamma[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$v_\gamma$")
plt.title(r"$v_\gamma$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/v_gamma.png"
)

# plot v_cdm and v_b together
plt.figure(figsize=(8, 6))
plt.plot(x, v_cdm[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$", color="C0")
plt.plot(x, v_b[0], color="C0", linestyle="--")
plt.plot(x, v_cdm[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$", color="C1")
plt.plot(x, v_b[1], color="C1", linestyle="--")
plt.plot(x, v_cdm[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$", color="C2")
plt.plot(x, v_b[2], color="C2", linestyle="--")
plt.xlabel(r"$x$")
plt.ylabel(r"$v$")
plt.yscale("log")
plt.title(r"$v_{CDM}$ and $v_{b}$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/v_cdm_v_b.png"
)

# plot v_nu
plt.figure(figsize=(8, 6))
plt.plot(x, v_nu[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, v_nu[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, v_nu[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.xlabel(r"$x$")
plt.ylabel(r"$v_\nu$")
plt.title(r"$v_\nu$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/v_nu.png"
)

# get photon and neutrino quadrupoles
Theta2 = np.zeros((3, npts))
Nu2 = np.zeros((3, npts))
for ki in range(len(k)):
    for xi in range(npts):
        Theta2[ki, xi] = pert.get_Theta(x[xi], k[ki], 2)
        Nu2[ki, xi] = pert.get_Nu(x[xi], k[ki], 2)

# old plot Theta2 and Nu2
plt.figure(figsize=(8, 6))
plt.plot(x, Theta2[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, Nu2[0], color="C0", linestyle="--")
plt.plot(x, Theta2[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, Nu2[1], color="C1", linestyle="--")
plt.plot(x, Theta2[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.plot(x, Nu2[2], color="C2", linestyle="--")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Theta_2$ and $\nu_2$")
plt.title(r"$\Theta_2$ and $\nu_2$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/Theta2_Nu2.png"
)
# plt.show()

# get gravitational potentials
Phi = np.zeros((3, npts))
Psi = np.zeros((3, npts))
for ki in range(len(k)):
    for xi in range(npts):
        Phi[ki, xi] = pert.get_Phi(x[xi], k[ki])
        Psi[ki, xi] = pert.get_Psi(x[xi], k[ki])

# plot Phi and Phi+Psi
plt.figure(figsize=(8, 6))
plt.plot(x, Phi[0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, Phi[0] + Psi[0], color="C0", linestyle="--")
plt.plot(x, Phi[1], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, Phi[1] + Psi[1], color="C1", linestyle="--")
plt.plot(x, Phi[2], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.plot(x, Phi[2] + Psi[2], color="C2", linestyle="--")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Phi$ and $\Phi+\Psi$")
plt.title(r"$\Phi$ and $\Phi+\Psi$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/Phi_Psi.png"
)

# get the polarization multipoles 0, 1 and 2
Theta_p = np.zeros((3, 3, npts))
for ell in range(3):
    for ki in range(len(k)):
        for xi in range(npts):
            Theta_p[ell, ki, xi] = pert.get_Theta_p(x[xi], k[ki], ell)

# plot Theta_p
plt.figure(figsize=(8, 6))
plt.plot(x, Theta_p[0, 0], label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(x, Theta_p[0, 1], color="C0", linestyle="--")
plt.plot(x, Theta_p[0, 2], color="C0", linestyle=":")
plt.plot(x, Theta_p[1, 0], label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(x, Theta_p[1, 1], color="C1", linestyle="--")
plt.plot(x, Theta_p[1, 2], color="C1", linestyle=":")
plt.plot(x, Theta_p[2, 0], label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.plot(x, Theta_p[2, 1], color="C2", linestyle="--")
plt.plot(x, Theta_p[2, 2], color="C2", linestyle=":")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Theta_p$")
plt.title(r"$\Theta_p$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/Theta_p.png"
)

# tests

# Theta[0] + Psi approx cos(k eta / sqrt(3))

x_end_rec = -7
x_rec = np.linspace(x_start, x_end_rec, npts)

# get Theta, Psi and eta before recombination
Theta_0 = np.zeros((3, npts))
Psi = np.zeros((3, npts))
eta = np.zeros((3, npts))

for ki in range(len(k)):
    for xi in range(npts):
        Theta_0[ki, xi] = pert.get_Theta(x_rec[xi], k[ki], 0)
        Psi = pert.get_Psi(x_rec[xi], k[ki])
        eta[ki, xi] = cosmo.eta_of_x(x_rec[xi])

# plot Theta_0 + Psi and cos(k eta / sqrt(3))
plt.figure(figsize=(8, 6))
plt.plot(x_rec, Theta_0[0] + Psi, label=r"$k=0.001 \mathrm{Mpc}^{-1}$")
plt.plot(
    x_rec,
    np.cos(k[0] * eta[0] / np.sqrt(3)),
    color="C0",
    linestyle="--",
)
plt.plot(x_rec, Theta_0[1] + Psi, label=r"$k=0.01 \mathrm{Mpc}^{-1}$")
plt.plot(
    x_rec,
    np.cos(k[1] * eta[1] / np.sqrt(3)),
    color="C1",
    linestyle="--",
)
plt.plot(x_rec, Theta_0[2] + Psi, label=r"$k=0.1 \mathrm{Mpc}^{-1}$")
plt.plot(
    x_rec,
    np.cos(k[2] * eta[2] / np.sqrt(3)),
    color="C2",
    linestyle="--",
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\Theta_0 + \Psi$ and $\cos(k \eta / \sqrt{3})$")
plt.title(r"$\Theta_0 + \Psi$ and $\cos(k \eta / \sqrt{3})$ vs $x$")
plt.legend()
plt.grid()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Perturbations/Theta_0_Psi_cos.png"
)
