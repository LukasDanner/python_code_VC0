import numpy as np
import sys

from scipy import special
from scipy.interpolate import interp1d

from helpers_classical import lookup_in_dict
from helpers_classical import fourier_transform
from helpers_classical import inverse_fourier_transform

class Circuit_VCOclSp:
    '''

    '''

    def __init__(self,
                 w0 = 1.0, Qcoil = 10.0, alpha_od = 1.0, n = 1.3, r_eff = 0.0005, # VCO parameters 
                 wL = 1.0, T1 = 5000, T2=5000, K = 0.0,                          # spin dynamics
                 Vn = 0.0,
                 **params):

        # VCO parameters

        # damping rate and quality factor:
        # Qcoil = w0/gamma_t = w0 L / R
        self.Qcoil = Qcoil

        print(self.Qcoil)

        self.gamma_t = 1.0 / self.Qcoil
    
        # eigenfrequency scale
        self.w0 = w0

        # transistor slope
        self.n = n

        # overdrive parameter alpha_od = G_mo / (2 G_t)
        self.alpha_od = alpha_od

        # effective geometry of coil used 
        self.r_eff = r_eff

        # SPIN parameters
        
        # Larmor frequency, can not be zero due to normalization to M0
        self.wL = wL

        # spin relaxation times
        self.T1 = T1
        self.T2 = T2

        # couling constant between oscillator and spins
        self.K = K  

        # initial values 
        self.initial = np.array([
            lookup_in_dict(params, "y0", 0.0),
            lookup_in_dict(params, "y1", 0.0),
            lookup_in_dict(params, "Mx", 0.0),
            lookup_in_dict(params, "My", 0.0),
            lookup_in_dict(params, "Mz", 1.0),
            0.0 # initial current is always zero
        ])

        # voltage noise 
        self.Vn = Vn

        self.include_noise = True if self.Vn != 0.0 else False
        
        print("build rhs of system")
        self.rhs = self.build_rhs()



    def build_rhs(self):

        res = np.zeros(len(self.initial))

        def rhs(t, y):

            # vco dof
            res[0] = y[1]
            res[1] = - y[0]  - (1.0/self.Qcoil) * ((1 - self.alpha_od) + 3 * self.alpha_od /(8 * self.n**2) * y[0]**2) * y[1] - self.K * y[2]

            # spin dof 
            res[2] = - y[2] / self.T2 + self.wL * y[3]

            # coupling with vco
            res[3] = - y[3] / self.T2 - self.wL * y[2] + y[4] * self.r_eff * y[5]
            res[4] = - (y[4]-1.0) / self.T1 - y[3] * self.r_eff * y[5]

            # test line for cos-B field
            #res[3] = - y[3] / self.T2 - self.wL * y[2] + y[4] * self.r_eff * np.cos(self.w0 *t)
            #res[4] = - (y[4]-1.0) / self.T1 - y[3] * self.r_eff * np.cos(self.w0 *t)

            # coil current
            res[5] = self.Qcoil /(2 * self.alpha_od) * y[0]

            return res

        def rhs_switched(y, t):

            return rhs(t, y)

        return rhs_switched if self.include_noise else rhs
    
    def build_Gfunc(self):

        G_Noise = np.zeros([len(self.initial), len(self.initial)])
        G_Noise[1,1] = - self.Vn 

        def Gfunc(y, t):
            return G_Noise
        
        return Gfunc



