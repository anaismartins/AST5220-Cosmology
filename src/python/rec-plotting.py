import matplotlib.pyplot as plt
import numpy as np

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/recombination.txt"

data = np.loadtxt(file)

x = data[:, 0]
Xe = data[:, 1]
Xe_Saha = data[:, 2]
ne = data[:, 3]
tau = data[:, 4]
dtaudx = data[:, 5]
ddtauddx = data[:, 6]
g_tilde = data[:, 7]
dgdx = data[:, 8]
ddgddx = data[:, 9]

x_lss = -6.990913222773562
x_rec = -6.97011625341667
x_rec_saha = -7.140341934397215

plt.plot(x, Xe_Saha, label="Saha")
plt.plot(x, Xe, label="Peebles")
plt.axvline(x=x_lss, linestyle="--", label="Last Scattering Surface", color="black")
plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
plt.axvline(
    x=x_rec_saha, linestyle=":", label="Recombination according to Saha", color="black"
)
plt.xlabel(r"$x$")
plt.ylabel(r"$X_e$")
plt.yscale("log")
plt.ylim(1e-4, 10)
plt.tight_layout()
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/Xe.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, tau, label=r"$\tau$")
plt.plot(x, -dtaudx, label=r"$-\frac{d\tau}{dx}$")
plt.plot(x, ddtauddx, label=r"$\frac{d^2\tau}{dx^2}$")
# plot vertical line at x = x_lss
plt.axvline(x=x_lss, linestyle="--", label="Last Scattering Surface", color="black")
plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
plt.axvline(
    x=x_rec_saha, linestyle=":", label="Recombination according to Saha", color="black"
)
plt.xlabel(r"$x$")
plt.yscale("log")
plt.ylim(1e-8, 1e8)
plt.xlim(-12, 0)
plt.tight_layout()
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/tau.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, g_tilde)
plt.axvline(x=x_lss, linestyle="--", label="Last Scattering Surface", color="black")
plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
plt.axvline(
    x=x_rec_saha, linestyle=":", label="Recombination according to Saha", color="black"
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{g}$")
plt.legend()
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/g_tilde.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, dgdx)
plt.axvline(x=x_lss, linestyle="--", label="Last Scattering Surface", color="black")
plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
plt.axvline(
    x=x_rec_saha, linestyle=":", label="Recombination according to Saha", color="black"
)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\tilde{g}}{dx}$")
plt.legend()
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/dgdx.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, ddgddx)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d^2\tilde{g}}{dx^2}$")
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/ddgddx.pdf"
)
# plt.show()
plt.clf()
