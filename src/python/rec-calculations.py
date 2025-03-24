import numpy as np
from BackgroundCosmology import BackgroundCosmology
from RecombinationHistory import RecombinationHistory

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

cosmo = BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB)
cosmo.solve()

rec = RecombinationHistory(
    cosmo=cosmo,
    Yp=Yp,
    z_reion=z_reion,
    delta_z_reion=delta_z_reion,
    z_Hereion=z_Hereion,
    delta_z_Hereion=delta_z_Hereion,
    reion=False,
)
rec.solve()

x_start = np.log(1e-8)
x_end = 0

x = np.linspace(x_start, x_end, 1000000)

# last scattering surface
g_tilde = np.zeros_like(x)
for i in range(len(x)):
    g_tilde[i] = rec.g_tilde_of_x(x[i])

x_lss = x[np.argmax(g_tilde)]
z_lss = cosmo.get_z(x_lss)
t_lss = cosmo.get_cosmic_time(x_lss)

print(f"The last scattering surface is at x = {x_lss}, z = {z_lss}, t = {t_lss}")

# recombination
# find when Xe = 0.1
Xe = np.zeros_like(x)
for i in range(len(x)):
    Xe[i] = rec.Xe_of_x(x[i])
    if Xe[i] < 0.1:
        x_recombination = x[i]
        break

z_recombination = cosmo.get_z(x_recombination)
t_recombination = cosmo.get_cosmic_time(x_recombination)

print(
    f"Recombination happens at x = {x_recombination}, z = {z_recombination}, t = {t_recombination}"
)
