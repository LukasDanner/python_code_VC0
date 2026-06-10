import numpy as np
import sys

from scipy import special
from scipy.interpolate import interp1d

from helpers_classical import lookup_in_dict
from helpers_classical import fourier_transform
from helpers_classical import inverse_fourier_transform

class Circuit_VCOclSp_RF:
    '''

    '''

    def __init__(self,
                Qcoil = 10.0, alpha_od = 1.0, n = 1.3, r_eff = 0.0005, # VCO parameters 
                wL = 1.0, T1 = 5000, T2=5000, K = 0.0,                          # spin dynamics
                Vn = 0.0,
                **params):

        # RWA or not
        self.do_RWA = lookup_in_dict(params, "do_RWA", False) 

        # sum cutoff for the approximation of the Bx-field
        self.nMAX = lookup_in_dict(params, "nMAX", 0)

        # VCO parameters

        # damping rate and quality factor:
        # Qcoil = w0/gamma_t = w0 L / R
        self.Qcoil = Qcoil


        self.gamma_t = 1.0 / self.Qcoil

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


    def get_B_function(self):

        def B_function(t, y):

            fac = - self.Qcoil /(2 * self.alpha_od) * self.r_eff 
            res = 0

            for nn in range(self.nMAX):

                res += fac * (1.0j)**(nn+1)/self.Qcoil**nn * 


    def build_rhs(self):

        res = np.zeros(len(self.initial))

        Q = self.Qcoil

        fac = 3 * self.alpha_od /(8 * self.n**2)

        i = 1.0j

        T1tilde = Q * self.T1
        T2tilde = Q * self.T2

        Delta = Q * (self.wL -1)        

        if self.do_RWA:

        else:

            def rhs(t, y):

                # vco dof
                res[0] = y[1]

                res[1] = -2 * i * Q * y[1] - (y[1] + i * y[0]) * (1 - self.alpha_od + fac * np.abs(y[0])**2) 
                - fac * (y[0]**2 *np.conjugate(y[1]) + y[1] * np.abs(y[0])**2 + y[0]**2 * (y[1] + i * Q * y[0]) * np.exp(2*i*Q*t))


                # spin dof 
                res[2] = -y[2] * (1 / self.T2 + i * Delta ) + i * y[3] * Q * 

                res[3] = (1-y[3])/T1tilde - Q * np.imag(y[2] * )

                # B-field
                res[4] = - self.Qcoil /(2 * self.alpha_od) * self.r_eff * (i * y[0] - y[1]/Q)

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



