'''
Functions for creating the fortran gridsearch module.
'''

# Standard libraries
import os
import glob

# External libraries
import numpy as np
from numpy import f2py


def create():
    '''
    Gives the user the option to automatically create a Python C/API extension
    module for the gridsearch algorithm.

    If the module is created, the user should re-run SKHASH and the module
    should be automatically used.
    '''
    usr_input=input('The Fortran subroutine module does not exist, but you have \'use_fortran=True\'. Would you like to create this module? [y / n]\n')
    usr_input=usr_input.strip().lower()
    if (usr_input!='y') & (usr_input!='yes'):
        print('\nThe module will not be created. If you would like to create the module manually, the following may help.')
        print('\tStep 1. Navigate to the SKHASH/functions folder, e.g.:\n\t\tcd functions')
        print('\tStep 2. Run the following command:\n\t\tpython3 -m numpy.f2py -c gridsearch.f -m gridsearch')
        return 0
    else:
        print('\nCreating module. This may take a minute...')
        with open("functions/gridsearch.f") as sourcefile:
            sourcecode = sourcefile.read()
        result=f2py.compile(sourcecode, modulename='gridsearch')
        if result==0:
            so_filepaths=glob.glob('gridsearch.cpython-*-*.so')
            for so_filepath in so_filepaths:
                os.rename(so_filepath,os.path.join('functions',so_filepath))
            try:
                import functions.gridsearch as gridsearch
                print('\n\nThe fortran extension module has been successfully created!\n')
                return 1
            except Exception as e:
                raise ValueError('There was an issue loading the module. Error message:\n{}'.format(e))
        else:
            raise ValueError('The compiled extension module could not be created.')


def prep_subroutine(sr_azimuth,takeoff,p_pol,sp_amp,npick0,nmc):
    '''
    Gets the polarities, S/P ratios, azimuth, and takeoff information in the
    correct format to call the Fortan subroutine.
    '''
    p_azi_mc=np.zeros((npick0,500))
    p_the_mc=np.zeros((npick0,500))
    p_azi_mc[:sr_azimuth.shape[0],:nmc]=sr_azimuth
    p_the_mc[:sr_azimuth.shape[0],:nmc]=takeoff

    npsta=sr_azimuth.shape[0]

    f_sp_amp=np.zeros(npsta)
    f_sp_amp[:len(sp_amp)]=sp_amp
    f_sp_amp[np.isnan(f_sp_amp)]=-9

    f_p_pol=np.zeros(npsta,float)
    f_p_pol[:len(p_pol)]=p_pol
    p_qual=np.abs(f_p_pol).astype(float)
    f_p_pol[f_p_pol<0]=-1
    f_p_pol[f_p_pol>0]=1
    f_p_pol=f_p_pol.astype(int)

    return p_azi_mc,p_the_mc,f_sp_amp,f_p_pol,p_qual
