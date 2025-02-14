import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as const

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/cosmology.txt"

# Read the file
data = np.loadtxt(file)

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

# Plot
plt.plot(x, dHpdx/Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\mathcal{H}}{H'}$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/dHpdx_over_Hp.pdf")
plt.clf()

plt.plot(x, ddHpddx/Hp)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d^2\mathcal{H}}{dx^2}$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/ddHpddx_over_Hp.pdf")
plt.clf()

plt.plot(x, eta * Hp / const.c)
plt.xlabel(r"$x$")
plt.ylabel(r"$\eta \mathcal{H}/c$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/etaHp_over_c.pdf")
plt.clf()

plt.plot(x, Hp)
plt.yscale("log")
plt.xlabel(r"$x$")
plt.ylabel(r"$\mathcal{H}$")
plt.savefig("/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/Hp.pdf")
plt.clf()