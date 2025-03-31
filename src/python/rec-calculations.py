import astropy.units as u
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
# g_tilde = np.zeros_like(x)
# for i in range(len(x)):
#     g_tilde[i] = rec.g_tilde_of_x(x[i])

tau = np.zeros_like(x)
for i in range(len(x)):
    tau[i] = rec.tau_of_x(x[i])
    if tau[i] < 1:
        x_lss = x[i]
        i_lss = i
        break

# x_lss = x[np.argmax(g_tilde)]
z_lss = cosmo.get_z(x_lss)
t_lss = cosmo.get_cosmic_time(x_lss) * u.s

print(
    f"The last scattering surface is at x = {x_lss}, z = {z_lss}, t = {t_lss.to(u.Myr)}"
)

# recombination
# find when Xe = 0.1
Xe = np.zeros_like(x)
Xe_Saha = np.zeros_like(x)
for i in range(len(x)):
    Xe[i] = rec.Xe_of_x(x[i])
    if Xe[i] < 0.1:
        x_recombination = x[i]
        break

for i in range(len(x)):
    Xe_Saha[i] = rec.Xe_of_x_Saha(x[i])
    if Xe_Saha[i] < 0.1:
        x_recombination_Saha = x[i]
        break

z_recombination = cosmo.get_z(x_recombination)
t_recombination = cosmo.get_cosmic_time(x_recombination) * u.s

z_recombination_Saha = cosmo.get_z(x_recombination_Saha)
t_recombination_Saha = cosmo.get_cosmic_time(x_recombination_Saha) * u.s

print(
    f"Recombination happens at x = {x_recombination}, z = {z_recombination}, t = {t_recombination.to(u.Myr)}"
)
print(
    f"Recombination happens at x = {x_recombination_Saha}, z = {z_recombination_Saha}, t = {t_recombination_Saha.to(u.Myr)}according to Saha"
)

# freeze-out abundance
Xe_freeze_out = rec.Xe_of_x(0.0)

print(f"Freeze-out abundance is Xe = {Xe_freeze_out}")

# sound horizon
r_s = rec.s_of_x(x_lss) * u.m

print(f"Sound horizon at last scattering surface is r_s = {r_s.to(u.Mpc)}")

# time of reionization
for i in range(i_lss, len(x)):
    # print(f"{x[i]}: {rec.Xe_of_x(x[i])}")
    if rec.Xe_of_x(x[i]) > 0.5:
        x_reion = x[i]
        break

z_reion = cosmo.get_z(x_reion)
t_reion = cosmo.get_cosmic_time(x_reion) * u.s

print(f"Reionization happens at x = {x_reion}, z = {z_reion}, t = {t_reion.to(u.Myr)}")
