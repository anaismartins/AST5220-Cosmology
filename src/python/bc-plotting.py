import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from BackgroundCosmology import BackgroundCosmology

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/cosmology.txt"

# Read the file
data = np.loadtxt(file)

h = 0.67
OmegaB0 = 0.05
OmegaCDM0 = 0.267
OmegaK0 = 0
Neff = 3.046
TCMB = 2.7255 * u.K

cosmo = BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB.value)
cosmo.solve()

x_radiation_matter = -8.131921542438867
x_matter_lambda = -0.25581887374692047

m = 1
km = 1e3 * m
s = 1
Mpc = 3.08567758e22 * m
H0_over_h = 100 * km / s / Mpc
c = 2.99792458e8 * m / s

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
plt.vlines(
    x_radiation_matter,
    min(dHpdx / Hp),
    max(dHpdx / Hp),
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    x_matter_lambda,
    min(dHpdx / Hp),
    max(dHpdx / Hp),
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\mathcal{H}}{dx}\frac{1}{\mathcal{H}}$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/dHpdx_over_Hp.pdf"
)
plt.clf()

plt.plot(x, ddHpddx / Hp)
plt.vlines(
    x_radiation_matter,
    min(ddHpddx / Hp),
    max(ddHpddx / Hp),
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    x_matter_lambda,
    min(ddHpddx / Hp),
    max(ddHpddx / Hp),
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
# plt.plot(x, ddHpddx, label="dH2dx")
# plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d^2\mathcal{H}}{dx^2}\frac{1}{\mathcal{H}}$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/ddHpddx_over_Hp.pdf"
)
plt.clf()

plt.plot(x, eta * Hp / c)
plt.vlines(
    x_radiation_matter,
    min(eta * Hp / c),
    max(eta * Hp / c),
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    x_matter_lambda,
    min(eta * Hp / c),
    max(eta * Hp / c),
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta \mathcal{H}/c$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/etaHp_over_c.pdf"
)
plt.clf()

x_acceleration = -0.4794171164592207

plt.plot(x, Hp / H0_over_h)
plt.vlines(
    x_acceleration,
    min(Hp / H0_over_h),
    max(Hp / H0_over_h),
    label="Universe starts to accelerate",
    color="black",
    linestyle="--",
)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\mathcal{H} \left(\frac{100 km/s}{Mpc}\right)$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/Hp.pdf"
)
plt.clf()

plt.plot(x, (t * u.s).to(u.Gyr).value)
plt.vlines(
    -8.1319,
    min(t * u.s).to(u.Gyr).value,
    max(t * u.s).to(u.Gyr).value,
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    -0.25582,
    min(t * u.s).to(u.Gyr).value,
    max(t * u.s).to(u.Gyr).value,
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$t$ (Gyr)")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/t.pdf"
)
plt.clf()

plt.plot(x, eta / Mpc / c)
plt.vlines(
    x_radiation_matter,
    min(eta / Mpc / c),
    max(eta / Mpc / c),
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    x_matter_lambda,
    min(eta / Mpc / c),
    max(eta / Mpc / c),
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta/c \left(\text{Mpc}\frac{\text{s}}{\text{m}}\right)$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/eta_over_c.pdf"
)
plt.clf()

plt.plot(x, OmegaB + OmegaCDM, label=r"$\Omega_M$")
plt.plot(x, OmegaGamma + OmegaNu, label=r"$\Omega_R$")
plt.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Omega_i$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/Omegas.pdf"
)
plt.clf()

supernova_data = np.loadtxt(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/supernovadata.txt"
)
z_supernova = supernova_data[:, 0]
dL_supernova = supernova_data[:, 1]
error_supernova = supernova_data[:, 2]

# get dl for these zs from the cpp function
z = np.logspace(np.log10(min(z_supernova)), np.log10(max(z_supernova)), 1000)
dL = np.zeros_like(z)
for zi in range(len(z)):
    dL[zi] = cosmo.get_luminosity_distance_of_x(np.log(1 / (z[zi] + 1)))

plt.plot(z, (dL / Mpc / 1e3 / z), label="Model")
plt.errorbar(
    z_supernova,
    dL_supernova / z_supernova,
    yerr=error_supernova,
    fmt=".",
    label="Supernova data",
)
plt.xscale("log")
# plt.yscale("log")
plt.xlabel(r"$z$")
plt.ylabel(r"$d_L/z$ (Gpc)")
# plt.xlim(min(z_supernova), max(z_supernova))
# plt.ylim(min(dL_supernova) - 0.02, max(dL_supernova) + 2)
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/dL.pdf"
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
plt.plot([0, 1], [1, 0], "--", label=r"Flat Universe", color="black")
plt.xlabel(r"$\Omega_M$")
plt.ylabel(r"$\Omega_\Lambda$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/fitting.pdf"
)
plt.clf()

plt.hist(h_fit[sigma1] * 100, bins=100)
# plt.vlines(67, 0, 100, label="Planck 2018", color="black", linestyle="--")
plt.xlabel(r"$H_0$ (km/s/Mpc)")
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/histogram.pdf"
)
plt.clf()

# eras plot
x_start = np.log(1e-8)
x_end = 0

plt.plot(x, OmegaB + OmegaCDM, label=r"$\Omega_M$")
plt.plot(x, OmegaGamma + OmegaNu, label=r"$\Omega_R$")
plt.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
plt.vlines(
    x_radiation_matter,
    0,
    1,
    label="Radiation-matter equality",
    color="black",
    linestyle="--",
)
plt.vlines(
    x_matter_lambda,
    0,
    1,
    label="Matter-dark energy equality",
    color="black",
    linestyles="-.",
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\Omega_i$")
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/Omegas.pdf"
)


# # plot eta Hp / c
# plt.plot(x, eta * Hp / const.c)
# plt.xlabel(r"$x$")
# plt.ylabel(r"$\eta \mathcal{H}/c$")
# plt.xscale("log")
# plt.xlim(x_start, x_end)
# plt.yscale("log")
# # plt.ylim(1e-3, 1e2)
# plt.savefig(
#     "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/BackgroundCosmology/etaHp_over_c.pdf"
# )
