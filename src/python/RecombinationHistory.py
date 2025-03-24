from ctypes import c_bool, c_double, c_void_p, cdll

lib = cdll.LoadLibrary(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/bin/librec.so"
)


class RecombinationHistory(object):
    def __init__(
        self, cosmo, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, reion
    ):
        lib.RecombinationHistory_new.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
            c_bool,
        ]
        lib.RecombinationHistory_new.restype = c_void_p

        print("Yp = ", Yp)

        lib.RecombinationHistory_solve.argtypes = [c_void_p]
        lib.RecombinationHistory_solve.restype = c_void_p

        print("z_reion = ", z_reion)

        lib.RecombinationHistory_tau_of_x.argtypes = [c_void_p, c_double]
        lib.RecombinationHistory_tau_of_x.restype = c_double

        lib.RecombinationHistory_g_tilde_of_x.argtypes = [c_void_p, c_double]
        lib.RecombinationHistory_g_tilde_of_x.restype = c_double

        print("delta_z_reion = ", delta_z_reion)

        lib.RecombinationHistory_Xe_of_x.argtypes = [c_void_p, c_double]
        lib.RecombinationHistory_Xe_of_x.restype = c_double

        self.obj = lib.RecombinationHistory_new(
            cosmo.obj, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion, reion
        )

    def solve(self):
        return lib.RecombinationHistory_solve(self.obj)

    def tau_of_x(self, x):
        return lib.RecombinationHistory_tau_of_x(self.obj, x)

    def g_tilde_of_x(self, x):
        return lib.RecombinationHistory_g_tilde_of_x(self.obj, x)

    def Xe_of_x(self, x):
        return lib.RecombinationHistory_Xe_of_x(self.obj, x)
