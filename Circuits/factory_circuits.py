# all circuits
from Circuits.circuit_VCO import Circuit_VCO

from helpers_classical import lookup_in_dict

def factory_circuits(**params):

    name = lookup_in_dict(params, 'circuit_name', None)

    if name == 'VCO':

        obj = Circuit_VCO(**params)

    else:

        print('ciruit-factory here: No circuit class could be found!')
        obj = None

    return obj

