import numpy as np
from helpers_classical import lookup_in_dict


def factory_timegrid(**params):
    '''
    returns a numpy-array with discretized time steps as given
    '''

    # initial time
    t0 = lookup_in_dict(params, 't0', 0.0)

    # final time
    tf = lookup_in_dict(params, 'tf', 0.0)

    # time step
    dt = lookup_in_dict(params, 'dt', 0.0)

    # number of time-intervals
    Nint = int((tf - t0)/dt + dt/2.0)

    # array with discretized time
    times = np.linspace(t0, tf, Nint + 1)

    return times