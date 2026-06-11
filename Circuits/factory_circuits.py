# all circuits
from Circuits.circuit_VCO import Circuit_VCO
from Circuits.circuit_VCOclSp import Circuit_VCOclSp
from Circuits.circuit_VCOclSp_RF import Circuit_VCOclSp_RF

from helpers_classical import lookup_in_dict

def factory_circuits(**params):

    name = lookup_in_dict(params, 'circuit_name', None)

    if name == 'VCO':

        obj = Circuit_VCO(**params)

    elif name == "VCOclSp":

        obj = Circuit_VCOclSp(**params)


    elif name == "VCOclSp_RF":

        obj = Circuit_VCOclSp_RF(**params)

    else:

        print('ciruit-factory here: No circuit class could be found!')
        obj = None

    return obj

