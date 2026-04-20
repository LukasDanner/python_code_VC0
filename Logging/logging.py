import os

from os import listdir
from os.path import isfile, join

from helpers_classical import lookup_in_dict
from helpers_classical import flatten_nested_dicts
from helpers_classical import search_nested_dict

class Logging:

    def __init__(self,
                 param_obj,                   # the parameter object used for the code
                 output_path='Hamiltonian/',  # output folder, seen relative to the python-project folder
                 **params):
        '''
        this object keeps track of storing results and automatically creates folders where files shall be stored.

        seen from the working directory of the main-file...
        -the root folder will be allocated at:
            output_path/circuit_name/folder_main/folder_run/
        -the sub-folders of single loop-runs will be allocated in the root-folder
         and have the name specified as parameter 'folder_vary' in the ini-file

        '''

        # dictionary of a single run (all flattened parameters)
        params_single_run = flatten_nested_dicts(param_obj.params, 0)

        # bools if we want to overwrite folders and file
        self.overwrite_folder_joined = lookup_in_dict(params_single_run, 'overwrite_folder_joined', True)
        self.overwrite_folder_single = lookup_in_dict(params_single_run, 'overwrite_folder_single', True)
        self.delete_existing_files = lookup_in_dict(params_single_run,  'delete_existing_files', True)

        # dictionary
        #       key    = name of the varying parameter
        #       value  = list of all the parameter values
        self.dict_vary = search_nested_dict(param_obj.params, param_obj.params['logging']['folder_vary'][0])

        # working directory of the code
        self.working_dir = os.getcwd()

        # output folder (relative to working directory)
        self.storage_dir_main = output_path

        # name of hamiltonian
        circuit_name = lookup_in_dict(params_single_run, 'circuit_name', None)

        # name of root-storage folder in output folder
        self.storage_dir_sub = lookup_in_dict(params_single_run, 'folder_main', None)

        # storage folder of this code run
        folder_run = lookup_in_dict(params_single_run, 'folder_run', 'DefaultRun')

        # complete path to root folder
        self.root = os.path.join(self.working_dir, self.storage_dir_main, circuit_name, self.storage_dir_sub, folder_run)

        # complete path to subfolder in root (-> is set later)
        self.subdir = None

        # create root directory
        self.root = self.handle_create_folder(self.root, 'root')

    def handle_create_folder(self, dir_name, dir_type):
        '''
        creates a folder (either 'root' or 'subfolder')
        while handling the case if we want to delete already exisiting subfolders and files before the loop-run

        returns the name of the created folder

        '''

        # bool if we want to overwrite folder
        overwrite_folder = self.overwrite_folder_joined if dir_type == 'root' else self.overwrite_folder_single

        # we can overwrite folders in this folder
        if overwrite_folder:

            # create folder
            self.create_dir(dir_name)

        # we do not want to overwrite folders in the possibly existing folder
        else:

            # if directory exists
            # then find a new name
            if os.path.exists(dir_name):

                # a test name for the new directory
                test = dir_name

                # loop over possible new names
                for ii in range(1, 1000):

                    # new folder name
                    test = dir_name + '_' + str(ii)

                    # if new name does not exist, exit loop
                    if not os.path.exists(test):

                        # found a new directory name
                        break

                # create the renamed directory
                self.create_dir(test)

                # set the directory name which was created
                dir_name = test

            # if directory does not exist, then just create it
            else:

                # create folder if it does not exist yet
                self.create_dir(dir_name)

        # have to delete all text files of old run in root folder
        if dir_type == 'root' and self.overwrite_folder_joined and self.delete_existing_files:

            print('logging here: we have deleted all files in root folder (but not in its sub-folders)')
            print('              we will overwrite subfolders if same name exists')

            # get all file names
            onlyfiles = [f for f in listdir(dir_name) if isfile(join(dir_name, f))]

            # remove these files
            for file in onlyfiles:
                os.remove(join(dir_name, file))

        return dir_name

    def create_dir(self, dirname):
        '''
        creates an unexisting folder (with full path dirname)
        '''

        if not os.path.exists(dirname):

            print('logging here: create the following directory')
            print('  '+dirname)
            os.makedirs(dirname)
    def get_alternative_keyname(self, key_name, circuit_name):
        '''
        return alternative key names for keys which have a unnecessary long name
        '''

        if key_name == "dim_H1":
            
            return "d1"
        
        if key_name == "dim_H2":
            
            return "d2"
        
        if key_name == "delta1":
            
            return "D1"
      
        if key_name == "delta2":
            
            return "D2"
        
        if key_name == "zpf1":
            
            return "z1"

        if key_name == "zpf2":
            
            return "z2"     
        
        if key_name == "xratio":
            
            return "x"

        return None

    def notify_new_run(self, N_now, separator='=', format_str='{:.3f}', circuit_name=""):
        '''
        if this function is called, then we create a new sub-folder,
        in which results/files for the current loop-run will be stored
        '''

        # the name of the varied parameter
        dummy1 = list(self.dict_vary.keys())[0]

        # if there exists an alterative key-name, then take it
        dummy2 = self.get_alternative_keyname(dummy1, circuit_name)

        # the parameter-name of the new sub-folder
        subdir_name = dummy1 if dummy2 is None else dummy2

        # the current value (of the current loop-run) of the varied parameter
        subdir_value = self.dict_vary[dummy1][N_now]

        # the full name of the sub-folder
        if type(subdir_value) is str:
            subdir_str = subdir_name + separator + subdir_value
        else:
            subdir_str = subdir_name + separator + format_str.format(subdir_value)

        # the full path of the sub-folder
        self.subdir = os.path.join(self.root, subdir_str)

        # create the sub-folder
        self.subdir = self.handle_create_folder(self.subdir, 'subdir')