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

        # if the simulation uses protons or electrons
        self.simulate_protons = lookup_in_dict(params, "simulate_protons", 0)

        # VCO parameters

        # damping rate and quality factor:
        # Qcoil = w0/gamma_t = w0 L / R
        self.Qcoil = Qcoil

        # transistor slope
        self.n = n

        # overdrive parameter alpha_od = G_mo / (2 G_t)
        self.alpha_od = alpha_od

        # effective geometry of coil used 
        self.r_eff = r_eff

        # SPIN parameters
        
        # Larmor frequency, can not be zero due to normalization to M0
        self.wL = wL

        self.Delta = self.Qcoil * (self.wL -1) 

        # spin relaxation times
        self.T1 = T1 * self.Qcoil
        self.T2 = T2 * self.Qcoil

        # couling constant between oscillator and spins
        self.K = K  

        # initial values 
        self.initial = np.array([
            (lookup_in_dict(params, "y0", 0.0) -1.0j *  lookup_in_dict(params, "y1", 0.0))/2.0,
            0.0 + 0.0j,
            lookup_in_dict(params, "Mx", 0.0) + 1.0j * lookup_in_dict(params, "My", 0.0),
            lookup_in_dict(params, "Mz", 1.0),
            0.0 + 0.0j # initial B-field is always zero
        ])

        # voltage noise 
        self.Vn = Vn

        self.include_noise = True if self.Vn != 0.0 else False
        
        print("build rhs of system")
        self.rhs = self.build_rhs()


    def build_rhs(self):

        res = np.zeros(len(self.initial))

        Q = self.Qcoil

        fac = 3 * self.alpha_od /(8 * self.n**2)

        i = 1.0j

        pm = 1.0 if self.simulate_protons else -1.0


        if self.do_RWA:

            ### RWA
            def rhs(t, y):

                mxy_der = -y[2] * (1 / self.T2 + i * self.Delta ) + i * y[3] * Q * y[4]

                aval = np.conjugate(y[0]) if self.simulate_protons else y[0]
                mxy_val = np.conjugate(mxy_der) + i * Q * np.conjugate(y[2]) if self.simulate_protons else mxy_der + i * Q * y[2] 

                # vco dof
                res[0] = - 0.5 * y[0] * (1 - self.alpha_od + fac * np.abs(y[0])**2) + i * self.K * mxy_val / 4.0 * mxy_val

                res[1] = 0

                # spin dof 
                res[2] = mxy_der

                res[3] = (1-y[3])/self.T1 - Q * np.imag(y[2] * np.conjugate(y[4]))

                # B-field
                res[4] = pm * Q * (Q/(2 * self.alpha_od) * aval + i * y[4]) 

                return res
            
        ### FULL EOMs
        else:

            def rhs(t, y):

                mxy_der = -y[2] * (1 / self.T2 + i * self.Delta ) + i * y[3] * Q * y[4]

                aval = np.conjugate(y[0]) if self.simulate_protons else y[0]
                mxy_val = Q * np.conjugate(mxy_der) + i * Q**2 * np.conjugate(y[2]) if self.simulate_protons else Q* mxy_der + i * Q**2 * y[2] 

                # vco dof
                res[0] = y[1]

                res[1] = -2 * i * Q * y[1] - (y[1] + i * y[0]) * (1 - self.alpha_od + fac * np.abs(y[0])**2) 
                - fac * (y[0]**2 *np.conjugate(y[1]) + y[1] * np.abs(y[0])**2 + y[0]**2 * (y[1] + i * Q * y[0]) * np.exp(2*i*Q*t))
                + self.K * mxy_val / 2.0

                # spin dof 
                res[2] = mxy_der

                res[3] = (1-y[3])/self.T1 - Q * np.imag(y[2] * np.conjugate(y[4]))

                # B-field
                res[4] = pm * Q * (Q/(2 * self.alpha_od) * aval + i * y[4]) 

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



