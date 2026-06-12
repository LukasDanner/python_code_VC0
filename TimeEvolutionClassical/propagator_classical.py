import sys
import numpy as np
from scipy.integrate import ode, complex_ode, solve_ivp
from sdeint import itoEuler


from helpers_classical import lookup_in_dict
from helpers_classical import log_to_file


class PropagatorClassical:

    def __init__(self, times, circuit,
                 print_propagation_status=True,
                 integrator='dop853', method='bdf',
                 atol=1.49012e-8, rtol=1.49012e-8, nsteps=1000,
                 include_noise = False, t_every=1,
                 **params):

        self.t_every = int(t_every + 0.1)
        self.include_noise = include_noise

        # object of the hamiltonian class
        self.times = times

        # dimensions of the state, given as an array
        self.circuit = circuit

        # results-vector
        # will store all results at all discrete times
        self.res = None

        # details for propagation
        self.print_propagation_status = print_propagation_status

        # if only a steady state shall be computed
        # (only makes sense for time-independent hamiltonians)
        self.do_only_steady_state = lookup_in_dict(params, 'do_only_steady_state', False)

        # name of the integrator and the used method
        self.integrator = integrator
        self.method = method

        # absolute and relative integration tolerance
        self.atol = float(atol)
        self.rtol = float(rtol)

        # number of max. internal integration steps of the python-integrator
        self.nsteps = nsteps

        self.solver_options = {"nsteps": self.nsteps,
                               "atol": self.atol,
                               "rtol": self.rtol,
                               "method": self.method,
                               "progress_bar": "text"
                               }


    def allocate_results(self, COMPLEX=False):
        '''
        allocate memory-space for propagation results

        - self.mat shall be an array which has the length of the number of discrete time steps
        - at each array index kk, self.mat[kk] shall be a matrix with full Hilbert space dimension
        '''

        LEN = int(len(self.times) / self.t_every) + 1 if self.t_every > 1 else len(self.times)

        # if we have not allocated a matrix yet
        if self.res is None or len(self.res) < LEN:

            if COMPLEX:


                # allocate space for this matrix
                self.res = np.empty([LEN, len(self.circuit.initial)], dtype=np.complex_)

            else:

                # allocate space for this matrix
                self.res = np.empty([LEN, len(self.circuit.initial)])

    def print_progress(self, t_cur, tf):
        '''
        quickly prints the status of the propagation
        '''

        # if print-statements shall be printed
        if self.print_propagation_status:
            print('     propagation is ' + '{:0.3f}'.format(100 * t_cur / tf) + '% completed')


    def propagate_ode(self, COMPLEX=False):
        '''
        - propagates an eom-system
        - stores the results in an array

        see also:
        https://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.integrate.ode.html
        '''

        # allocate solution matrix
        self.allocate_results(COMPLEX=COMPLEX)

        sol = solve_ivp(self.circuit.rhs, (self.times[0], self.times[-1]), self.circuit.initial, t_eval=self.times,
                        #rtol=self.rtol, atol=self.atol,
                        method=self.method
                        )

        self.res = sol.y.T

        print('status report of solver: ')    
        print(sol.success)
        print(sol.message)
        #print(sol.nfev)
        #print(sol.njev)
        #print(sol.nlu)
        #print(sol.t[-1], self.times[-1])

        '''

        # get an integrator
        if self.integrator == "zvode":

            r = ode(self.circuit.rhs).set_integrator(self.integrator,
                                                     method=self.method,
                                                     rtol=self.rtol,
                                                     atol=self.atol,
                                                     nsteps=self.nsteps)

        else:

            r = complex_ode(self.circuit.rhs).set_integrator(self.integrator,
                                                             method=self.method,
                                                             rtol=self.rtol,
                                                             atol=self.atol,
                                                             nsteps=self.nsteps)



        # set initial conditions
        r.set_initial_value(self.circuit.initial, self.times[0])

        # store the initial conditions for the solutions
        self.res[0] = self.circuit.initial


        # set the counter for the current time step
        step = 1

        # propagate
        while r.successful() and step < len(self.times):

            # integrate one step
            res = r.integrate(self.times[step])
            
            if step % self.t_every == 0:
                
                # store the result
                self.res[step] = res
                #self.circuit.store_sols(step, res)

            # print the progress
            if (100*step/(len(self.times) - 1)) % 10 == 0:

                self.print_progress(self.times[step], self.times[-1])

            # increase time-step
            step += 1

        # warn if propagation was unsuccessful
        if step < len(self.times):

            print('ATTENTION: propagation unsuccesful')

        '''

    def propagate_noise(self):

        Gfunc = self.circuit.build_Gfunc()

        self.allocate_results()

        self.res = itoEuler(self.circuit.rhs, Gfunc, self.circuit.initial, self.times)



    def find_steady_state(self):

        print("steady-state search not implemented")

    def propagate(self, **params):
        '''
        performs the time-evolution

        stores all results in
        '''

        COMPLEX = lookup_in_dict(params, "COMPLEX", False)

        # only search for steady-state
        if self.do_only_steady_state:

            self.find_steady_state()

        # make full time evolution
        else:


            if self.include_noise:

                print('propagating with sdeint-itoEuler: ')

                self.propagate_noise()

            else:

                if self.integrator == 'odeintw':

                    print('ATTENTION: NOT IMPLEMENTED - INTEGRATION FAILED')
                    #self.propagate_odeintw()

                else:
                    
                    #print('propagating with integrator: ' + self.integrator)
                    self.propagate_ode(COMPLEX=COMPLEX)

