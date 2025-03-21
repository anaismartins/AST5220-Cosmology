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

plt.plot(x, Xe_Saha, label="Saha")
plt.plot(x, Xe, label="Peebles", linestyle="--")
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
plt.plot(x, -dtaudx, label=r"$-\frac{d\tau}{dx}$", linestyle="--")
plt.plot(x, ddtauddx, label=r"$\frac{d^2\tau}{dx^2}$", linestyle=":")
plt.xlabel(r"$x$")
plt.yscale("log")
plt.ylim(1e-8, 1e8)
plt.tight_layout()
plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/tau.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, g_tilde)
plt.xlabel(r"$x$")
plt.ylabel(r"$\tilde{g}$")
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/g_tilde.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, dgdx)
plt.xlabel(r"$x$")
plt.ylabel(r"$\frac{d\tilde{g}}{dx}$")
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
