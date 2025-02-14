import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/cosmology.txt"

# Read the file
data = np.loadtxt(file)

m = 1
km = 1e3 * m
s = 1
Mpc = 3.08567758e22 * m
H0_over_h = 100 * km/s/Mpc

# Extract the columns
x = data[:,0]
eta = data[:,1]
Hp = data[:,2]
dHpdx = data[:,3]
ddHpddx = data[:,4]
OmegaB = data[:,5]
OmegaCDM = data[:,6]
OmegaLambda = data[:,7]
OmegaGamma = data[:,8]
OmegaNu = data[:,9]
OmegaK = data[:,10]
t = data[:,11]
z = data[:,12]
dL = data[:,13]

# Plot
plt.plot(x, dHpdx/Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\mathcal{H}}{dx}\frac{1}{\mathcal{H}}$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/dHpdx_over_Hp.pdf")
plt.clf()

plt.plot(x, ddHpddx/Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d^2\mathcal{H}}{dx^2}\frac{1}{\mathcal{H}}$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/ddHpddx_over_Hp.pdf")
plt.clf()

plt.plot(x, eta * Hp / const.c)
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta \mathcal{H}/c$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/etaHp_over_c.pdf")
plt.clf()

plt.plot(x, Hp/H0_over_h)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\mathcal{H} \left(\frac{100 km/s}{Mpc}\right)$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Hp.pdf")
plt.clf()

plt.plot(x, t)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$t$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/t.pdf")
plt.clf()

plt.plot(x, eta/const.c)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta/c$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/eta_over_c.pdf")
plt.clf()

plt.plot(x, OmegaB + OmegaCDM, label=r"$\Omega_M$")
plt.plot(x, OmegaGamma + OmegaNu, label=r"$\Omega_R$")
plt.plot(x, OmegaLambda, label=r"$\Omega_\Lambda$")
plt.xlabel(r"$x$")
plt.ylabel(r"$\Omega_i$")
plt.legend()
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Omegas.pdf")
plt.clf()

supernova_data = np.loadtxt("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/data/supernovadata.txt")
z_supernova = supernova_data[:,0]
dL_supernova = supernova_data[:,1]
error_supernova = supernova_data[:,2]

plt.plot(z, dL)
plt.errorbar(z_supernova, dL_supernova, yerr=error_supernova, fmt='o')
plt.yscale("log")
plt.xlabel(r"$z$")
plt.ylabel(r"$d_L$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/dL.pdf")
plt.clf()