import configparser
import math
import sys

import numpy as np

from helpers_classical import str_to_array
from helpers_classical import change_datatype_array_elements
from helpers_classical import max_length_array_in_dict
from helpers_classical import arr_to_arr_of_given_length

class Parameters:

    def __init__(self):

        # file name of input file -> will be set later
        self.ini_fname = None

        # number of total runs -> will be set later
        self.N_runs = 1

        # the run number for this loop run
        self.N_now = 0

        # this is a dictionary of dictionaries, where...
        #   ...the keys are the section names of the ini file
        #   ...the value at each key is a dictionary itself, where...
        #       ...the keys are the parameter names
        #       ...the values at each key are arrays of parameter values for this parameter name
        self.params = dict()


        # this is a dictionary of dictionaries, where...
        #   ...the keys are the section names of the ini file
        #   ...the value at each key is a dictionary itself, where...
        #       ...the keys are the parameter names
        #       ...the value at each key is the parameter value for this loop-run
        self.params_now = dict()

        # this is a flattened version of the dictionary self.params_now
        self.params_now_flat = dict()

        # get a list of section names of the ini file
        self.ini_keys = None

    def setup_all_parameters(self, fname):
        '''
        - reads the ini-file
        - stores all parameters given
        '''

        # set the file name of the input file
        self.ini_fname = fname

        # read the input file
        self.params, self.ini_keys = self.read_input(fname)

        # find the number of loop runs
        self.N_runs = self.find_Nruns(self.params)

        # adjust the length of the parameter lists
        self.params = self.adjust_length_of_arrays(self.params, self.N_runs)

    def read_input(self, fname):
        '''
        - reads the ini-file with name 'fname' into a nested dictionary
        - modifies the data-type of all given parameters

        - return the nested dictionary and an array containing all section-names of the ini-file
        '''

        # set up a parser
        config = configparser.ConfigParser(inline_comment_prefixes="#")

        # all parameters from config file are strings per default
        config.optionxform = str

        # read the complete config file
        config.read(fname)

        # the nested dictionary, as described as in the __init__ function
        param_dict = dict()

        # go through all sections from the ini-file
        #   - make arrays from input strings
        #   - change/modify data types
        #   - store all input parameters
        for section_name in config.sections():

            # a temporary dictionary at of the section
            tmp_dict = dict(config[section_name])

            # go through the temporary dictionary
            for key in tmp_dict:

                # make array from a string, where different inputs are separated by the separator
                tmp_dict[key] = str_to_array(tmp_dict[key], separator=',')

                # all input parameters of the 'timegrid'-section are doubles
                if section_name == 'timegrid':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'approximations'-section are integers
                if section_name == 'approximations':

                    # change the data-type of the array to integer-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'int')

                # all input parameters of the 'internal'-section are doubles
                if section_name == 'internal':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'noise'-section are doubles
                if section_name == 'noise':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'inputmode'-section are doubles
                if section_name == 'inputmode':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'outputmode'-section are doubles
                if section_name == 'outputmode':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')


                # all input parameters of the 'vectors'-section are itself arrays
                if section_name == 'vectors':

                    # change the data-type of the array to arrays
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'array')

                # all input parameters of the 'conditional'-section are bools
                if section_name == 'conditional':

                    # change the data-type of the array to bools
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'bool')

                # all input parameters of the 'default'-section are doubles
                if section_name == 'default':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'initial'-section are doubles
                if section_name == 'initial':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

                # all input parameters of the 'analysis'-section are doubles
                if section_name == 'analysis':

                    # change the data-type of the array to double-variables
                    tmp_dict[key] = change_datatype_array_elements(tmp_dict[key], 'np.double')

            # key: name of the section
            # value: the modified temporary dictionary
            param_dict[section_name] = tmp_dict

        return param_dict, config.sections()

    def find_Nruns(self, param_dict):

        # this will be the number of loops in the main-file
        N_runs = 0

        # go through all sections from the ini-file
        for section_name in param_dict:

            # get the length of the arrays in the sub-dictionaries
            test = max_length_array_in_dict(param_dict[section_name])

            # if this section contains arrays of larger length, modify N_runs
            if test > N_runs:

                N_runs = test

        # return the number of loop runs
        return N_runs

    def adjust_length_of_arrays(self, param_dict, N_runs):
        '''
        input:
            - a nested dictionary 'param_dict'
            - the number of loop-runs 'N_runs'

        output:
            - the modified 'param_dict', where all values of the sub-dictionaries
              will be arrays with length N_runs
        '''

        # go through all sections from the ini-file (= keys in the dictionary param_dict)
        # the values at each key is again a dictionary
        for section_name in param_dict:

            # go through all parameter names (=keys of the sub-dictionaries)
            # the values are arrays with different lengths (depending on the input parameters, separated by a comma)
            for key in param_dict[section_name]:

                # eventually make a bigger array for the value at the keys of the sub-dictionaries
                # if the array does not have length N_runs, then multiply the first array value N_runs times
                # else: leave the array like that
                param_dict[section_name][key] = arr_to_arr_of_given_length(param_dict[section_name][key], N_runs)

        return param_dict

    def setup_current_parameters(self, N_now):
        '''
        input:
            - the number of the current loop-run

        this function sets up the currently used parameters.

        '''

        # modify up loop run number
        self.N_now = N_now

        # go through all sections of the ini-file (=keys of the dictionary self.params)
        # the values at the keys are dictionaries itself
        for section_name in self.params:

            # create an empty dictionary at each section_name
            self.params_now[section_name] = dict()

            # go through the inner dictionary (key are all parameter-names of this section)
            # the values at the keys are arrays with the parameter values
            for key in self.params[section_name]:

                # the currently used parameter values are at position N_now in the arrays
                # we create the key and value of this
                self.params_now[section_name][key] = self.params[section_name][key][N_now]

                # we create a flattened dictionary
                # forgetting about the section names
                # and storing all keys+values (of the current loop-run) of the inner dictionary
                self.params_now_flat[key] = self.params[section_name][key][N_now]
