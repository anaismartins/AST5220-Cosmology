from BackgroundCosmology import BackgroundCosmology
from RecombinationHistory import RecombinationHistory
from Perturbations import Perturbations
from PowerSpectrum import PowerSpectrum
import numpy as np
import matplotlib.pyplot as plt

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

rec = RecombinationHistory(cosmo, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, reion)
rec.solve()

pert = Perturbations(cosmo, rec)
pert.solve()

A_s = 2.1e-9
n_s = 0.965
kpivot_mpc = 0.05

ps = PowerSpectrum(cosmo, rec, pert, A_s, n_s, kpivot_mpc)
ps.solve()

# plot temperature source function
# plt.figure(figsize=(10, 6))
# xs = np.linspace(-10, 0, 100)
# H0 = cosmo.get_H_of_x(0)
# ST = np.zeros(len(xs))
# for xi in range(len(xs)):
#     ST[xi] = pert.get_source_T(xs[xi], 340 * H0)
# plt.plot(xs, ST, label=r'$S_T$')
# plt.xlabel(r'$x$')
# plt.ylabel(r'$S_T$')
# plt.title(r'$S_T$')
# plt.axhline(0, color='black', lw=0.5)
# plt.axvline(0, color='black', lw=0.5)
# plt.grid()
# plt.legend()
# plt.show()

# get eta0 
eta0 = cosmo.eta_of_x(0)
H0 = cosmo.get_H_of_x(0)

m = 1.0
s = 1.0
Mpc = 3.08567758e22 * m
k_min = 0.00005 / Mpc
deltak = 2 * np.pi / eta0 / 6
k_max = 3000. / eta0
n_k = int((k_max - k_min) / deltak)

c = 2.99792458e8 * m / s

ks = np.linspace(k_min, k_max, n_k)

ells = np.array([2, 3, 4, 5, 6, 7, 8, 10, 12, 15,
        20, 25, 30, 40, 50, 60, 70, 80, 90, 100,
        120, 140, 160, 180, 200, 225, 250, 275, 300, 350,
        400, 450, 500, 550, 600, 650, 700, 750, 800, 850,
        900, 950, 1000, 1050, 1100, 1150, 1200, 1250, 1300, 1350,
        1400, 1450, 1500, 1550, 1600, 1650, 1700, 1750, 1800, 1850,
        1900, 1950, 2000])

thetaT = np.zeros((len(ells), len(ks)))

magic_ells = [6, 100, 200, 500, 1000]
magic_ellis = []

plt.figure(figsize=(10, 6))

for elli in range(len(ells)):
    ell = ells[elli]
    for ki in range(len(ks)):
        # print(f'ell = {ells[elli]}, k = {ks[ki] * Mpc}') 
        thetaT[elli, ki] = ps.get_thetaT_ell_of_k(ell, ks[ki])
        # print(f'ell = {ells[elli]}, k = {ks[ki] * Mpc}, thetaT = {thetaT[elli, ki]}')

    if ell in magic_ells:
        plt.plot(ks * c / H0, thetaT[elli, :], label=f'$\ell$ = {ell}')
        magic_ellis.append(elli)

plt.legend()
# plt.xlim(0, 500)
# plt.ylim(-0.015, 0.015)
plt.ylabel(r'$\Theta_\ell$')
plt.xlabel(r'$k c / H_0$')
plt.show()

# plot theta^2 / k
plt.figure(figsize=(10, 6))
for elli in magic_ellis:
    plt.plot(c * ks / H0, thetaT[elli, :]**2 * H0 / (ks * c) , label=f'$\ell$ = {ells[elli]}')

plt.xlabel(r'$k c / H_0$')
plt.ylabel(r'$\Theta_\ell^2 H_0/ (kc)$')
plt.legend()
plt.show()

# TODO: to debug, set S = g and get theta_l = (theta_0 + psi)_ini * f_l(k eta_0) (should be just the spherical bessel functions and l(l+1) C_l is like a constant)

ells = np.arange(2, 2000)
# plot power spectrum
plt.figure(figsize=(10, 6))
cell_TT = np.zeros(len(ells))
cell_plot = np.zeros(len(ells))
for elli in range(len(ells)):
    cell_TT[elli] = ps.get_cell_TT(ells[elli])
    cell_plot[elli] = ells[elli] * (ells[elli] + 1) / (2 * np.pi) * cell_TT[elli]

plt.plot(ells, cell_plot * 10e6**2)
plt.xlabel(r'$\ell$')
plt.ylabel(r'$\frac{\ell(\ell+1)C_\ell^{TT}}{2\pi}$ ($\mu K^2$)')
plt.xscale('log')
plt.yscale('log')
plt.show()

# plot matter power spectrum
plt.figure(figsize=(10, 6))
x = 0
matter_power_spectrum = np.zeros(len(ks))
for ki in range(len(ks)):
    matter_power_spectrum[ki] = ps.get_matter_power_spectrum(x, ks[ki])

matter_power_spectrum = matter_power_spectrum / (Mpc / h)**3
plt.plot(ks * Mpc / h, matter_power_spectrum)
plt.xlabel(r'$k$ (h/Mpc)')
plt.ylabel(r'$P(k)$')
plt.xscale('log')
plt.yscale('log')
plt.show()