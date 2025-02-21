import astropy.constants as const
import matplotlib.pyplot as plt
import numpy as np

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/cosmology.txt"

# Read the file
data = np.loadtxt(file)

m = 1
km = 1e3 * m
s = 1
Mpc = 3.08567758e22 * m
H0_over_h = 100 * km / s / Mpc

# Extract the columns
x = data[:, 0]
eta = data[:, 1]
Hp = data[:, 2]
dHpdx = data[:, 3]
ddHpddx = data[:, 4]
OmegaB = data[:, 5]
OmegaCDM = data[:, 6]
OmegaLambda = data[:, 7]
OmegaGamma = data[:, 8]
OmegaNu = data[:, 9]
OmegaK = data[:, 10]
t = data[:, 11]
z = data[:, 12]
dL = data[:, 13]

# Plot
plt.plot(x, dHpdx / Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\mathcal{H}}{dx}\frac{1}{\mathcal{H}}$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/dHpdx_over_Hp.pdf"
)
plt.clf()

plt.plot(x, ddHpddx / Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d^2\mathcal{H}}{dx^2}\frac{1}{\mathcal{H}}$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/ddHpddx_over_Hp.pdf"
)
plt.clf()

plt.plot(x, eta * Hp / const.c)
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta \mathcal{H}/c$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/etaHp_over_c.pdf"
)
plt.clf()

plt.plot(x, Hp / H0_over_h)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\mathcal{H} \left(\frac{100 km/s}{Mpc}\right)$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Hp.pdf"
)
plt.clf()

plt.plot(x, t)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$t$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/t.pdf"
)
plt.clf()

plt.plot(x, eta / const.c)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta/c$")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/eta_over_c.pdf"
)
plt.clf()

plt.plot(x, OmegaB + OmegaCDM, label=r"$\Omega_M$")
plt.plot(x, OmegaGamma + OmegaNu, label=r"$\Omega_R$")
plt.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Omega_i$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Omegas.pdf"
)
plt.clf()

supernova_data = np.loadtxt(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/supernovadata.txt"
)
z_supernova = supernova_data[:, 0]
dL_supernova = supernova_data[:, 1]
error_supernova = supernova_data[:, 2]

plt.plot(z, dL / Mpc / 1e3, label="Model")
plt.errorbar(
    z_supernova, dL_supernova, yerr=error_supernova, fmt=".", label="Supernova data"
)
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$z$")
plt.ylabel(r"$d_L$ (Gpc)")
plt.xlim(min(z_supernova) - 0.002, max(z_supernova) + 0.2)
plt.ylim(min(dL_supernova) - 0.02, max(dL_supernova) + 2)
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/dL.pdf"
)
plt.clf()

fit_data = np.loadtxt(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/results_supernovafitting.txt",
    skiprows=200,
)
chi2_fit = fit_data[:, 0]
h_fit = fit_data[:, 1]
OmegaM_fit = fit_data[:, 2]
OmegaK_fit = fit_data[:, 3]

OmegaLambda_fit = (
    1 - OmegaM_fit - OmegaK_fit
)  # assuming OmegaGamma = 0 since we are today

chi2_norm = chi2_fit - min(chi2_fit)
sigma1 = chi2_norm < 3.53

plt.plot(
    OmegaM_fit[sigma1], OmegaLambda_fit[sigma1], ".", label=r"$1\sigma$ constraint"
)
# line where OmegaM = OmegaLambda
plt.plot([0, 1], [1, 0], "--", label=r"$\Omega_M = \Omega_\Lambda$", color="black")
plt.xlabel(r"$\Omega_M$")
plt.ylabel(r"$\Omega_\Lambda$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/fitting.pdf"
)
plt.clf()

plt.hist(h_fit[sigma1] * 100, bins=100)
plt.xlabel(r"$H_0$ (km/s/Mpc)")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/histogram.pdf"
)
plt.clf()
