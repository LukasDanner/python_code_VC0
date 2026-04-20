import numpy as np
import sys

from scipy import special
from scipy.interpolate import interp1d

from helpers_classical import lookup_in_dict
from helpers_classical import fourier_transform
from helpers_classical import inverse_fourier_transform

class Circuit_VCO:
    '''

    '''

    def __init__(self,
                 w0 = 1.0, Qfac = 10.0, alpha_od = 1.0, n = 1.3, # system parameters 
                 Vn = 0.0,
                 eta =0.0, 
                 **params):


        # damping rate and quality factor:
        # Qcoil = w0/gamma_t = w0 L / R
        self.Qcoil = Qfac
        self.gamma_t = 1.0 / self.Qcoil
    
        # eigenfrequency scale
        self.w0 = w0

        # transistor slope
        self.n = n

        # overdrive parameter alpha_od = G_mo / (2 G_t)
        self.alpha_od = alpha_od

        # initial values 
        self.initial = np.array([
            lookup_in_dict(params, "y0", 0.0),
            lookup_in_dict(params, "y1", 0.0),
        ])

        # voltage noise 
        self.Vn = Vn

        # filling factor of ESR material inside coil
        self.eta = eta 

        self.include_noise = True if self.Vn != 0.0 else False
        
        print("build rhs of system")
        self.rhs = self.build_rhs()



    def build_rhs(self):

        res = np.zeros(len(self.initial))

        def rhs(t, y):

            res[0] = y[1]

            res[1] = - y[0]  - (1.0/self.Qcoil) * ((1 - self.alpha_od) + 3 * self.alpha_od /(8 * self.n) * y[0]**2) * y[1]

            return res

        def rhs_switched(y, t):

            return rhs(t, y)

        return rhs_switched if self.include_noise else rhs
    
    def build_Gfunc(self):

        G_Noise = np.zeros([len(self.initial), len(self.initial)])
        G_Noise[-1,-1] = - self.Vn 

        def Gfunc(y, t):
            return G_Noise
        
        return Gfunc



