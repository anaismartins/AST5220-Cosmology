import astropy.constants as const
import astropy.units as u
import numpy as np
from BackgroundCosmology import BackgroundCosmology

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

# radiation-matter equality
print("Radiation-matter equality")

a = OmegaGamma0 / OmegaR0
x = np.log(a)
print(f"x = {x}")

z = 1 / np.exp(x) - 1
print(f"z = {z}")


t_radiation_matter = cosmo.get_cosmic_time(x) * u.s
print(f"t = {t_radiation_matter.to(u.Gyr)}")
