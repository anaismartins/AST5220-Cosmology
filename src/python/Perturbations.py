from ctypes import c_double, c_void_p, cdll

lib = cdll.LoadLibrary(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/bin/libcmb.so"
)


class Perturbations(object):
    def __init__(self, cosmo, rec):
        lib.Perturbations_new.argtypes = [c_void_p, c_void_p]
        lib.Perturbations_new.restype = c_void_p

        lib.Perturbations_solve.argtypes = [c_void_p]
        lib.Perturbations_solve.restype = c_void_p

        lib.Perturbations_get_delta_cdm.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_delta_cdm.restype = c_double

        lib.Perturbations_get_delta_b.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_delta_b.restype = c_double

        lib.Perturbations_get_Theta.argtypes = [c_void_p, c_double, c_double, c_double]
        lib.Perturbations_get_Theta.restype = c_double

        lib.Perturbations_get_Nu.argtypes = [c_void_p, c_double, c_double, c_double]
        lib.Perturbations_get_Nu.restype = c_double

        lib.Perturbations_get_v_cdm.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_v_cdm.restype = c_double

        lib.Perturbations_get_v_b.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_v_b.restype = c_double

        lib.Perturbations_get_Phi.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_Phi.restype = c_double

        lib.Perturbations_get_Psi.argtypes = [c_void_p, c_double, c_double]
        lib.Perturbations_get_Psi.restype = c_double

        lib.Perturbations_get_Theta_p.argtypes = [
            c_void_p,
            c_double,
            c_double,
            c_double,
        ]
        lib.Perturbations_get_Theta_p.restype = c_double

        self.obj = lib.Perturbations_new(cosmo.obj, rec.obj)

    def solve(self):
        return lib.Perturbations_solve(self.obj)

    def get_delta_cdm(self, x, k):
        return lib.Perturbations_get_delta_cdm(self.obj, x, k)

    def get_delta_b(self, x, k):
        return lib.Perturbations_get_delta_b(self.obj, x, k)

    def get_Theta(self, x, k, ell):
        return lib.Perturbations_get_Theta(self.obj, x, k, ell)

    def get_Nu(self, x, k, ell):
        return lib.Perturbations_get_Nu(self.obj, x, k, ell)

    def get_v_cdm(self, x, k):
        return lib.Perturbations_get_v_cdm(self.obj, x, k)

    def get_v_b(self, x, k):
        return lib.Perturbations_get_v_b(self.obj, x, k)

    def get_Phi(self, x, k):
        return lib.Perturbations_get_Phi(self.obj, x, k)

    def get_Psi(self, x, k):
        return lib.Perturbations_get_Psi(self.obj, x, k)

    def get_Theta_p(self, x, k, ell):
        return lib.Perturbations_get_Theta_p(self.obj, x, k, ell)
