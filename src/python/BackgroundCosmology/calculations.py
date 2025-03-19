import astropy.constants as const
import astropy.units as u
import numpy as np
from BackgroundCosmology import BackgroundCosmology

m = 1
Mpc = 3.08567758e22 * m

h = 0.67
OmegaB0 = 0.05
OmegaCDM0 = 0.267
OmegaK0 = 0
Neff = 3.046
TCMB = 2.7255 * u.K

cosmo = BackgroundCosmology(h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB.value)
cosmo.solve()

OmegaGamma0 = cosmo.get_OmegaGamma(x=0)
OmegaNu0 = cosmo.get_OmegaNu(x=0)
OmegaR0 = OmegaGamma0 + OmegaNu0

OmegaM0 = OmegaB0 + OmegaCDM0

OmegaLambda0 = 1 - OmegaM0 - OmegaR0 - OmegaK0

# radiation-matter equality
print("Radiation-matter equality")

a = OmegaR0 / OmegaM0
x = np.log(a)
print(f"x = {x}")

z = 1 / np.exp(x) - 1
print(f"z = {z}")

t_radiation_matter = cosmo.get_cosmic_time(x) * u.s
print(f"t = {t_radiation_matter.to(u.Gyr)}")

# matter-dark energy equality
print("Matter-dark energy equality")

a = np.power(OmegaM0 / OmegaLambda0, 1 / 3)
x = np.log(a)
print(f"x = {x}")

z = 1 / np.exp(x) - 1
print(f"z = {z}")

t_matter_lambda = cosmo.get_cosmic_time(x) * u.s
print(f"t = {t_matter_lambda.to(u.Gyr)}")

# when the universe starts to accelerate
print("Universe starts to accelerate")

# \ddot a = 0
# H Hp' = 0

x_start = np.log(1e-8)
x_end = 0

x_vec = np.linspace(x_start, x_end, 1000)
H_vec = np.zeros_like(x_vec)
dHpdx_vec = np.zeros_like(x_vec)

for i, x in enumerate(x_vec):
    H_vec[i] = cosmo.get_H_of_x(x)
    dHpdx_vec[i] = cosmo.get_dHpdx_of_x(x)

x_acceleration = x_vec[np.argmin(np.abs(dHpdx_vec))]
print(f"x = {x_acceleration}")

z = 1 / np.exp(x_acceleration) - 1
print(f"z = {z}")

t_acceleration = cosmo.get_cosmic_time(x_acceleration) * u.s
print(f"t = {t_acceleration.to(u.Gyr)}")

# conformal time today
print("Conformal time today")
x_today = 0
eta_today = cosmo.eta_of_x(x_today) / Mpc / 1e3
print(f"eta = {eta_today} Gpc")
