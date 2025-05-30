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

x_lss = -6.994726307500645
x_rec = -6.97011625341667
x_rec_saha = -7.140341934397215
x_reion = -2.2014024950700097

# plt.plot(x, Xe_Saha, label="Saha")
# plt.plot(x, Xe, label="Saha + Peebles")
plt.plot(x, Xe)

# plt.axvline(x=x_lss, linestyle="--", label="LSS", color="black")
# plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
# plt.axvline(x=x_rec_saha, linestyle=":", label="Recombination (Saha)", color="black")

plt.hlines(
    y=0.1,
    xmin=min(x),
    xmax=(x_lss),
    color="C2",
)
plt.text(min(x) + 1, 0.125, "Plasma", fontsize=10, color="C2")
plt.hlines(y=0.1, xmin=x_lss, xmax=x_reion, color="C3")
plt.text(x_lss + 0.3, 0.125, "Neutral\nhydrogen\nand helium", fontsize=10, color="C3")
plt.hlines(y=0.1, xmin=x_reion, xmax=0, color="C4")
plt.text(
    x_reion + 0.15, 0.125, "Neutral\nand ionised\nspecies", fontsize=10, color="C4"
)


plt.xlabel(r"$x$")
plt.ylabel(r"$X_e$")
# plt.yscale("log")
# plt.ylim(1e-4, 10)
# plt.xlim(-10, -4)
plt.xlim(min(x), max(x) + 1)
plt.tight_layout()
# plt.legend()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/Xe.pdf"
)
# plt.show()
plt.clf()

plt.plot(x, tau, label=r"$\tau$")
plt.plot(x, -dtaudx, label=r"$-\frac{d\tau}{dx}$")
plt.plot(x, ddtauddx, label=r"$\frac{d^2\tau}{dx^2}$")
# plot vertical line at x = x_lss
plt.axvline(x=x_lss, linestyle="--", label="Last Scattering", color="black")
# plt.axvline(x=x_rec, linestyle="-.", label="Recombination", color="black")
# plt.axvline(
#     x=x_rec_saha, linestyle=":", label="Recombination according to Saha", color="black"
# )
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

plt.plot(x, g_tilde / max(abs(g_tilde)), label=r"$\tilde{g}$")
plt.plot(x, dgdx / max(abs(dgdx)), label=r"$\frac{d\tilde{g}}{dx}$")
plt.plot(x, ddgddx / max(abs(ddgddx)), label=r"$\frac{d^2\tilde{g}}{dx^2}$")
plt.xlabel(r"$x$")
plt.legend()
plt.xlim(-8, -1.5)
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/gtilde.pdf"
)
plt.show()
plt.clf()
