import matplotlib.pyplot as plt
import numpy as np

file = "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/data/recombination.txt"

data = np.loadtxt(file)

x = data[:, 0]
Xe = data[:, 1]

plt.plot(x, Xe)
plt.xlabel(r"$x$")
plt.ylabel(r"$X_e$")
# plt.yscale("log")
plt.tight_layout()
plt.savefig(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/output/plots/RecombinationHistory/Xe.pdf"
)
plt.clf()
