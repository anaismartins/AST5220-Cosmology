from ctypes import c_double, c_void_p, cdll

lib = cdll.LoadLibrary(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/bin/libcmb.so"
)


class RecombinationHistory(object):
    def __init__(self, Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion):
        lib.RecombinationHistory_new.argtypes = [
            c_double,
            c_double,
            c_double,
            c_double,
            c_double,
        ]
        lib.RecombinationHistory_new.restype = c_void_p

        lib.RecombinationHistory_solve.argtypes = [c_void_p]
        lib.RecombinationHistory_solve.restype = c_double

        lib.RecombinationHistory_tau_of_x.argtypes = [c_void_p, c_double]
        lib.RecombinationHistory_tau_of_x.restype = c_double

        self.obj = lib.RecombinationHistory_new(
            Yp, z_reion, delta_z_reion, z_Hereion, delta_z_Hereion
        )

    def solve(self):
        return lib.RecombinationHistory_solve(self.obj)

    def tau_of_x(self, x):
        return lib.RecombinationHistory_tau_of_x(self.obj, x)
