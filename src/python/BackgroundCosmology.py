from ctypes import c_double, c_void_p, cdll

lib = cdll.LoadLibrary(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/bin/libcmb.so"
)


class BackgroundCosmology(object):
    def __init__(self, h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB):
        lib.BackgroundCosmology_new.argtypes = [
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
        ]
        lib.BackgroundCosmology_new.restype = c_void_p

        lib.BackgroundCosmology_solve.argtypes = [c_void_p]
        lib.BackgroundCosmology_solve.restype = c_double

        lib.BackgroundCosmology_get_cosmic_time.argtypes = [c_void_p, c_double]
        lib.BackgroundCosmology_get_cosmic_time.restype = c_double

        lib.BackgroundCosmology_get_OmegaGamma.argtypes = [c_void_p, c_double]
        lib.BackgroundCosmology_get_OmegaGamma.restype = c_double

        lib.BackgroundCosmology_get_OmegaNu.argtypes = [c_void_p, c_double]
        lib.BackgroundCosmology_get_OmegaNu.restype = c_double

        self.obj = lib.BackgroundCosmology_new(
            h, OmegaB0, OmegaCDM0, OmegaK0, Neff, TCMB
        )

    def solve(self):
        return lib.BackgroundCosmology_solve(self.obj)

    def get_cosmic_time(self, x):
        return lib.BackgroundCosmology_get_cosmic_time(self.obj, x)

    def get_OmegaGamma(self, x):
        return lib.BackgroundCosmology_get_OmegaGamma(self.obj, x)

    def get_OmegaNu(self, x):
        return lib.BackgroundCosmology_get_OmegaNu(self.obj, x)
