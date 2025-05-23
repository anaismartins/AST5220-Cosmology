from ctypes import c_double, c_void_p, cdll

lib = cdll.LoadLibrary(
    "/mn/stornext/u3/aimartin/d5/cosmologyii/AST5220-Cosmology/bin/libcmb.so"
)

class PowerSpectrum(object):
    def __init__(self, cosmo, rec, pert, A_s, n_s, kpivot_mpc):
        lib.PowerSpectrum_new.argtypes = [c_void_p, c_void_p, c_void_p, c_double, c_double, c_double]
        lib.PowerSpectrum_new.restype = c_void_p

        lib.PowerSpectrum_solve.argtypes = [c_void_p]
        lib.PowerSpectrum_solve.restype = c_void_p

        lib.PowerSpectrum_get_thetaT_ell_of_k.argtypes = [c_void_p, c_double, c_double]
        lib.PowerSpectrum_get_thetaT_ell_of_k.restype = c_double

        lib.PowerSpectrum_get_cell_TT.argtypes = [c_void_p, c_double]
        lib.PowerSpectrum_get_cell_TT.restype = c_double

        lib.PowerSpectrum_get_matter_power_spectrum.argtypes = [c_void_p, c_double, c_double]
        lib.PowerSpectrum_get_matter_power_spectrum.restype = c_double

        lib.PowerSpectrum_get_cell_EE.argtypes = [c_void_p, c_double]
        lib.PowerSpectrum_get_cell_EE.restype = c_double

        lib.PowerSpectrum_get_cell_TE.argtypes = [c_void_p, c_double]
        lib.PowerSpectrum_get_cell_TE.restype = c_double

        self.obj = lib.PowerSpectrum_new(cosmo.obj, rec.obj, pert.obj, A_s, n_s, kpivot_mpc)

    def solve(self):
        return lib.PowerSpectrum_solve(self.obj)
    
    def get_thetaT_ell_of_k(self, x, k):
        return lib.PowerSpectrum_get_thetaT_ell_of_k(self.obj, x, k)
    
    def get_cell_TT(self, ell):
        return lib.PowerSpectrum_get_cell_TT(self.obj, ell)
    
    def get_matter_power_spectrum(self, x, k):
        return lib.PowerSpectrum_get_matter_power_spectrum(self.obj, x, k)
    
    def get_cell_EE(self, ell):
        return lib.PowerSpectrum_get_cell_EE(self.obj, ell)
    
    def get_cell_TE(self, ell):
        return lib.PowerSpectrum_get_cell_TE(self.obj, ell)