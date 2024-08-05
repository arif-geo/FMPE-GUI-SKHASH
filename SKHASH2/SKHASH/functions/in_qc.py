'''
Functions used to provide quality control of user inputs
'''

# Standard libraries
import os
import multiprocessing

# External libraries
import numpy as np
import pandas as pd


def check_input_params(p_dict,qual_criteria_dict):
    '''
    Performs various checks to ensure that the user-defined parameters are appropriate
    '''
    # Ensures input format is lowercase and is not empty
    if (p_dict['input_format']==''):
        p_dict['input_format']='skhash'
    else:
        p_dict['input_format']=p_dict['input_format'].lower()

    for format_file in ['input_format_stfile','input_format_fpfile','input_format_impfile','input_format_conpfile','input_format_dlpfile']:
        p_dict[format_file]=p_dict[format_file].lower()

    if (p_dict['input_format'][:-1]=='hash') & (p_dict['catfile']!=''):
        raise ValueError('When using the {} input format, no catalog file (catfile) should be provided.'.format(p_dict['input_format']))

    if (p_dict['input_format'][:-1]=='hash') & (p_dict['fpfile']==''):
        raise ValueError('When using the {} input format, a p-polarity file (fpfile) should be provided.'.format(p_dict['input_format']))

    if (p_dict['input_format']=='hash1') | (p_dict['input_format']=='hash5'):
        if len(p_dict['vmodel_paths'])>0:
            raise ValueError('When using the {} input format, no velocity model(s) (vmodel_paths) should be provided.'.format(p_dict['input_format']))

    if (p_dict['input_format']=='hash1'): # hash1 only matches on the station and channel codes
        p_dict['require_network_match']=False
        p_dict['require_station_match']=True
        p_dict['require_location_match']=False
        p_dict['require_channel_match']=True
    elif p_dict['input_format'] in ['hash2','hash3','hash4','hash5']: # hash2-5 only matches on the network, station, and channel codes
        p_dict['require_network_match']=True
        p_dict['require_station_match']=True
        p_dict['require_location_match']=False
        p_dict['require_channel_match']=True

    default_input_format=p_dict['input_format']

    # Ensures the stfile format is an accepted format
    if p_dict['stfile']:
        possible_format=['skhash','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_stfile']=='':
            p_dict['input_format_stfile']=default_input_format
        if not(p_dict['input_format_stfile'] in possible_format):
            raise ValueError('input_format_stfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the fpfile format is an accepted format
    if p_dict['fpfile']:
        possible_format=['skhash','ncsn','hypoinverse','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_fpfile']=='':
            p_dict['input_format_fpfile']=default_input_format
        if not(p_dict['input_format_fpfile'] in possible_format):
            raise ValueError('input_format_fpfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the impfile format is an accepted format
    if p_dict['impfile']:
        possible_format=['skhash','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_impfile']=='':
            p_dict['input_format_impfile']=default_input_format
        if not(p_dict['input_format_impfile'] in possible_format):
            raise ValueError('input_format_impfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the conpfile format is an accepted format
    if p_dict['conpfile']:
        possible_format=['skhash','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_conpfile']=='':
            p_dict['input_format_conpfile']=default_input_format
        else:
            p_dict['input_format_conpfile']=p_dict['input_format_conpfile'].lower()
        if not(p_dict['input_format_conpfile'] in possible_format):
            raise ValueError('input_format_conpfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the dlpfile format is an accepted format
    if p_dict['dlpfile']:
        possible_format=['skhash','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_dlpfile']=='':
            p_dict['input_format_dlpfile']=default_input_format
        else:
            p_dict['input_format_dlpfile']=p_dict['input_format_dlpfile'].lower()
        if not(p_dict['input_format_dlpfile'] in possible_format):
            raise ValueError('input_format_dlpfile must be one of the following:\n\t{}'.format(possible_format))
        
    # Ensures the plfile format is an accepted format
    if p_dict['plfile']:
        possible_format=['skhash','hash1','hash2','hash3','hash4','hash5']
        if p_dict['input_format_plfile']=='':
            p_dict['input_format_plfile']=default_input_format
        else:
            p_dict['input_format_plfile']=p_dict['input_format_plfile'].lower()
        if not(p_dict['input_format_plfile'] in possible_format):
            raise ValueError('input_format_plfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the ampfile format is an accepted format
    if p_dict['ampfile']:
        possible_format=['skhash','hash3']
        if (p_dict['input_format_ampfile']==''):
            p_dict['input_format_ampfile']=default_input_format
        else:
            p_dict['input_format_ampfile']=p_dict['input_format_ampfile'].lower()
        if not(p_dict['input_format_ampfile'] in possible_format):
            raise ValueError('input_format_ampfile must be one of the following:\n\t{}'.format(possible_format))

    # Ensures the relative ampfile format is an accepted format
    if p_dict['relampfile']:
        possible_format=['skhash','hash3']
        if (p_dict['input_format_relampfile']==''):
            p_dict['input_format_relampfile']=default_input_format
        else:
            p_dict['input_format_relampfile']=p_dict['input_format_relampfile'].lower()
        if not(p_dict['input_format_relampfile'] in possible_format):
            raise ValueError('input_format_relampfile must be one of the following:\n\t{}'.format(possible_format))

    # Determines what station information should be used to associate picks with metadata.
    p_dict['merge_on']=[]
    if p_dict['require_network_match']:
        p_dict['merge_on'].append('network')
    if p_dict['require_station_match']:
        p_dict['merge_on'].append('station')
    if p_dict['require_location_match']:
        p_dict['merge_on'].append('location')
    if p_dict['require_channel_match']:
        p_dict['merge_on'].append('channel')

    # Ensures that require_location_match=False if doing HASH format as they do not consider location codes
    if (p_dict['input_format'][:-1]=='hash') & (p_dict['require_location_match']):
        print('*WARNING: Metadata cannot be associated using location codes (require_location_match=True) when using the {} input_format. Setting require_location_match=False.'.format(p_dict['input_format']))
        p_dict['require_location_match']=False

    if p_dict['input_format']=='skhash':
        if not(p_dict['require_network_match'] | p_dict['require_station_match'] | p_dict['require_location_match'] | p_dict['require_channel_match']):
            raise ValueError('When using the SKHASH input format, at least one of the following metadata matches must be True: require_network_match, require_station_match, require_location_match, require_channel_match')

    if (not(p_dict['fpfile'])) & (not(p_dict['impfile'])) & (not(p_dict['conpfile'])) & (not(p_dict['dlpfile'])):
        print('*WARNING: No P polarity files provided (fpfile, impfile, conpfile, dlpfile), so there will be no polarity agreement output file (outfile_pol_agree: {}).'.format(p_dict['outfile_pol_agree']))
        p_dict['outfile_pol_agree']=''
    if (not(p_dict['ampfile'])) & (not(p_dict['relampfile'])) & (p_dict['corfile']!=''):
        print('*WARNING: No traditional or consensus S/P amplitude file provided (ampfile or relampfile), so the S/P correction file (corfile: {}) will be ignored.'.format(p_dict['corfile']))
        p_dict['corfile']=''

    # Ensures all provided input files exist
    if p_dict['catfile']:
        if not(os.path.isfile(p_dict['catfile'])):
            raise ValueError('Catalog file (catfile) does not exist:{}'.format(p_dict['catfile']))
    if p_dict['stfile']:
        if not(os.path.isfile(p_dict['stfile'])):
            raise ValueError('Station file (stfile) does not exist:{}'.format(p_dict['stfile']))
    if p_dict['plfile']:
        if not(os.path.isfile(p_dict['plfile'])):
            raise ValueError('Polarity reversal file (plfile) does not exist:{}'.format(p_dict['plfile']))
    if p_dict['corfile']:
        if not(os.path.isfile(p_dict['corfile'])):
            raise ValueError('Station correction file (corfile) does not exist:{}'.format(p_dict['corfile']))
    if p_dict['ampfile']:
        if not(os.path.isfile(p_dict['ampfile'])):
            raise ValueError('S/P amplitude file (ampfile) does not exist:{}'.format(p_dict['ampfile']))
    if p_dict['relampfile']:
        if not(os.path.isfile(p_dict['relampfile'])):
            raise ValueError('Relative S/P amplitude file (relampfile) does not exist:{}'.format(p_dict['relampfile']))
    if p_dict['fpfile']:
        if not(os.path.isfile(p_dict['fpfile'])):
            raise ValueError('P-polarity file (fpfile) does not exist:{}'.format(p_dict['fpfile']))
    if p_dict['impfile']:
        if not(os.path.isfile(p_dict['impfile'])):
            raise ValueError('Imputed P-polarity file (impfile) does not exist:{}'.format(p_dict['impfile']))
    if p_dict['conpfile']:
        if not(os.path.isfile(p_dict['conpfile'])):
            raise ValueError('Consensus P-polarity file (conpfile) does not exist:{}'.format(p_dict['conpfile']))
    # Ensures the velocity model paths also exist:
    for tmp_vmodel_path in p_dict['vmodel_paths']:
        if not(os.path.isfile(tmp_vmodel_path)):
            raise ValueError('Velocity model (vmodel_paths) path does not exist: {}'.format(tmp_vmodel_path))

    # Ensures the output folder for plots exists
    if p_dict['outfolder_plots']:
        if os.path.exists(p_dict['outfolder_plots']):
            if not(os.path.isdir(p_dict['outfolder_plots'])):
                raise ValueError('Folder containg output plots (outfolder_plots) exists but is not a folder: {}'.format(p_dict['outfolder_plots']))

    # If plotting station names, ensures station names are being considered
    if p_dict['outfolder_plots']:
        if (p_dict['plot_station_names'] & (p_dict['require_station_match']==False)):
            raise ValueError('Station names can be added to plots ($plot_station_names==True) only if you are considering station names ($require_station_match==True)')

    # Ensures all of the input and output filepaths are unique
    filepath_vars=[
        'catfile',
        'stfile',
        'plfile',
        'corfile',
        'fpfile',
        'impfile',
        'conpfile',
        'ampfile',
        'relampfile',
        'simulpsfile',
        'outfile1',
        'outfile2',
        'outfile_pol_agree',
        'outfile_sp_agree',
        'outfile_pol_info',
        'outfolder_plots']
    tmp_dict = {key: p_dict[key] for key in filepath_vars if p_dict[key]!=''}
    if len(tmp_dict)!=len(set(tmp_dict.values())):
        rev_multidict = {}
        for key, value in tmp_dict.items():
            rev_multidict.setdefault(value, set()).add(key)
        tmp=[values for key, values in rev_multidict.items() if len(values) > 1]
        raise ValueError('Filepaths for the following variables are repeated: {}'.format(tmp[0]))

    # Ensures existing files are not overwritten by the output, if desired
    if (not(p_dict['overwrite_output_file'])):
        if p_dict['outfile1']:
            if os.path.exists(p_dict['outfile1']):
                raise ValueError('Preferred mechanism output file (outfile1={}) already exists. Either change the outfile1 path, remove this file, or set overwrite_output_file=True.'.format(p_dict['outfile1']))
        if p_dict['outfile2']:
            if os.path.exists(p_dict['outfile2']):
                raise ValueError('Acceptable mechanism output file (outfile2={}) already exists. Either change the outfile1 path, remove this file, or set overwrite_output_file=True.'.format(p_dict['outfile2']))
        if p_dict['outfile_pol_agree']:
            if os.path.exists(p_dict['outfile_pol_agree']):
                raise ValueError('Polarity agreement output file (outfile_pol_agree={}) already exists. Either change the outfile1 path, remove this file, or set overwrite_output_file=True.'.format(p_dict['outfile_pol_agree']))
        if p_dict['outfile_sp_agree']:
            if os.path.exists(p_dict['outfile_sp_agree']):
                raise ValueError('S/P difference output file (outfile_sp_agree={}) already exists. Either change the outfile1 path, remove this file, or set overwrite_output_file=True.'.format(p_dict['outfile_sp_agree']))
        if p_dict['outfile_pol_info']:
            if os.path.exists(p_dict['outfile_pol_info']):
                raise ValueError('Polarity info output file (outfile_pol_info={}) already exists. Either change the outfile1 path, remove this file, or set overwrite_output_file=True.'.format(p_dict['outfile_pol_info']))

    # Creates output directories if necessary
    if p_dict['outfile1']:
        folder_path=os.path.dirname(p_dict['outfile1'])
        if folder_path:
            os.makedirs(folder_path,exist_ok=True)
    if p_dict['outfile2']:
        folder_path=os.path.dirname(p_dict['outfile2'])
        if folder_path:
            os.makedirs(folder_path,exist_ok=True)
    if p_dict['outfile_pol_agree']:
        folder_path=os.path.dirname(p_dict['outfile_pol_agree'])
        if folder_path:
            os.makedirs(folder_path,exist_ok=True)
    if p_dict['outfile_sp_agree']:
        folder_path=os.path.dirname(p_dict['outfile_sp_agree'])
        if folder_path:
            os.makedirs(folder_path,exist_ok=True)
    if p_dict['outfile_pol_info']:
        folder_path=os.path.dirname(p_dict['outfile_pol_info'])
        if folder_path:
            os.makedirs(folder_path,exist_ok=True)
    if p_dict['outfolder_plots']:
        os.makedirs(p_dict['outfolder_plots'],exist_ok=True)

    if p_dict['min_amp']>1:
        raise ValueError('min_amp must be a value <=1')
    if (p_dict['max_agap']<0) | (p_dict['max_agap']>360):
        raise ValueError('The max azimuthal gap (max_agap) must be a value between 0 and 360')
    if (p_dict['max_pgap']<0) | (p_dict['max_pgap']>90):
        raise ValueError('The max takeoff angle gap (max_pgap) must be a value between 0 and 90')
    if p_dict['maxout']<1:
        raise ValueError('The max number of acceptable focal mechanisms (maxout) must be at least 1 (ideally larger!).')
    if p_dict['nmc']<1:
        raise ValueError('The number of trials (nmc) must be at least 1 (ideally larger!).')

    if p_dict['min_quality_report']:
        if not(p_dict['min_quality_report'] in qual_criteria_dict['qual_letter']):
            raise ValueError('The minimum mech quality (min_quality_report: {}) must be one of the quality codes:\n\t{}'.format(p_dict['min_quality_report'],qual_criteria_dict['qual_letter']))

    if p_dict['delmax']>p_dict['look_del'][1]:
        raise ValueError('The maximum source-receiver distance (delmax: {} km) is greater than the model distance range (look_del: {}-{} km)'.format(p_dict['delmax'],p_dict['look_del'][0],p_dict['look_del'][1]))

    # Makes sure the number of requested cores are available
    if p_dict['num_cpus']!=1:
        num_cpus_avail=multiprocessing.cpu_count()-1
        if (p_dict['num_cpus']==0): # Uses all available cores
            p_dict['num_cpus']=num_cpus_avail
        elif (p_dict['num_cpus']>num_cpus_avail): # If the num requested > available cores, use all available
            print('You requested {} cores, but only {} are available. Using all available cores instead.'.format(p_dict['num_cpus'],num_cpus_avail))
            p_dict['num_cpus']=num_cpus_avail

    if p_dict['outfolder_plots']:
        try:
            import matplotlib.pyplot as plt
        except:
            raise ImportError('matplotlib was unable to be loaded, but it is required to make plots. Either install or remove $outfolder_plots from your control file.')

    return p_dict


def check_cat_df(cat_df,pol_df,p_dict):
    '''
    Does quality control testing on the catalog dataframe (cat_df).
    '''
    if len(cat_df)==0:
        return pd.DataFrame()

    if not('vert_uncert_km' in cat_df.columns):
        cat_df['vert_uncert_km']=0
    if not('horz_uncert_km' in cat_df.columns):
        cat_df['horz_uncert_km']=0

    # Ensures the event_id column is string
    cat_df['event_id']=cat_df['event_id'].astype(str)

    # Rounds lat/lon to the specified number of decimal places
    if p_dict['epicenter_degree_precision']>0:
        cat_df['origin_lon']=cat_df['origin_lon'].round(p_dict['epicenter_degree_precision'])
        cat_df['origin_lat']=cat_df['origin_lat'].round(p_dict['epicenter_degree_precision'])

    # Ensures the vertical and horizontal uncertainties are positive
    if np.any(cat_df['vert_uncert_km']<0):
        prob_ind=cat_df['vert_uncert_km'].argmin()
        raise ValueError('Vertical uncertainties must be >= 0. Example problem: event {}, vert uncertainty: {} km'.
            format(cat_df.loc[prob_ind,'event_id'],cat_df.loc[prob_ind,'vert_uncert_km']))
    if np.any(cat_df['horz_uncert_km']<0):
        prob_ind=cat_df['horz_uncert_km'].argmin()
        raise ValueError('Horizontal uncertainties must be >= 0. Example problem: event {}, horz uncertainty: {} km'.
            format(cat_df.loc[prob_ind,'event_id'],cat_df.loc[prob_ind,'horz_uncert_km']))

    if p_dict['perturb_epicentral_location']:
        if (p_dict['default_vert_uncert_km']==0) & (cat_df['vert_uncert_km'].min()==0):
            print('*WARNING: You set default_vert_uncert_km=0 and there are events with 0 km uncertainties. The hypocentral perturbation will not work correctly.')
        else:
            novert_uncert_ind=np.where(cat_df['vert_uncert_km']==0)[0]
            if len(novert_uncert_ind)>0:
                print('{} events have a 0 km vertical uncertainty. Setting the uncertainty to {} km for these events'.
                    format(len(novert_uncert_ind),p_dict['default_vert_uncert_km']))
                cat_df.loc[novert_uncert_ind,'vert_uncert_km']=p_dict['default_vert_uncert_km']
        if (p_dict['default_horz_uncert_km']==0) & (cat_df['horz_uncert_km'].min()==0):
            print('*WARNING: You set default_horz_uncert_km=0 and there are events with 0 km uncertainties. The hypocentral perturbation will not work correctly.')
        else:
            nohorz_uncert_ind=np.where(cat_df['horz_uncert_km']==0)[0]
            if len(nohorz_uncert_ind)>0:
                print('{} events have a 0 km horizontal uncertainty. Setting the uncertainty to {} km for these events'.
                    format(len(nohorz_uncert_ind),p_dict['default_horz_uncert_km']))
                cat_df.loc[nohorz_uncert_ind,'horz_uncert_km']=p_dict['default_horz_uncert_km']

    # Selects only cataloged events that have a polarity measurement
    cat_df=cat_df.loc[cat_df['event_id'].isin(pol_df['event_id']),:].reset_index(drop=True)

    return cat_df


def check_pol_df(pol_df,cat_df,p_dict):
    '''
    Does quality control testing on the polarity dataframe (pol_df).
    '''
    # Ensures the event_id column is string
    pol_df['event_id']=pol_df['event_id'].astype(str)

    if not('p_polarity' in pol_df.columns):
        pol_df['p_polarity']=np.nan
    if not('sr_dist_km' in pol_df.columns):
        pol_df['sr_dist_km']=np.nan

    # Applies default location uncertainties if none were given by the user
    if not('horz_uncert_km' in pol_df.columns):
        pol_df['horz_uncert_km']=p_dict['default_horz_uncert_km']
    if not('vert_uncert_km' in pol_df.columns):
        pol_df['vert_uncert_km']=p_dict['default_vert_uncert_km']

    # Applies default location uncertainties if location uncertainties are <0
    pol_df.loc[ (pd.isnull(pol_df['horz_uncert_km'])) | (pol_df['horz_uncert_km']<0),'horz_uncert_km']=p_dict['default_horz_uncert_km']
    pol_df.loc[ (pd.isnull(pol_df['vert_uncert_km'])) | (pol_df['vert_uncert_km']<0),'vert_uncert_km']=p_dict['default_vert_uncert_km']

    # Checks to make sure all event_ids in the polarity file(s) are included in the catalog.
    if len(cat_df)>0:
        unique_pol_event_id=pol_df['event_id'].unique()
        id_found_flag=np.isin(unique_pol_event_id,cat_df['event_id'])
        if np.any(~id_found_flag):
            raise ValueError('{} events in polarity or SP files not found in catalog.\n\tMissing problematic IDs: {}'.
                format(np.sum(~id_found_flag),unique_pol_event_id[np.where(~id_found_flag)[0]]))

    # Ensures that there are not null values for the polarities/sp_ratios
    if 'sp_ratio' in pol_df:
        null_flag=pd.isnull(pol_df[['p_polarity','sp_ratio']]).all(axis=1)
        if null_flag.any():
            raise ValueError('Null values are present for both the polarity and S/P columns. Example problematic records:\n{}'.
                            format(pol_df.loc[null_flag,['event_id','sta_code','p_polarity','sp_ratio','source']].head(10).to_string()))
    else:
        null_flag=pd.isnull(pol_df['p_polarity'])
        if null_flag.any():
            raise ValueError('Null values are present for the polarity column. Example problematic records:\n{}'.
                            format(pol_df.loc[null_flag,['event_id','sta_code','p_polarity','source']].head(10).to_string()))

    # # Ensures there aren't null values in the polarity dataset
    # if pol_df.isnull().values.any():
    null_pol_df=pol_df.drop(pol_df.filter(['sr_dist_km','takeoff','azimuth','takeoff_uncertainty','azimuth_uncertainty','sp_ratio','p_polarity','station']),axis=1)
    if null_pol_df.isnull().values.any():
        null_index=null_pol_df[null_pol_df.isna().any(axis=1)].index
        if len(null_index)>5:
            null_index=null_index[:5]
        raise ValueError('There are null values in datasets. If you are using multiple polarity files, do the columns match?\n{}'.format(null_pol_df.loc[null_index,:].to_string()))

    # Ensures P polarity weights are in the appropriate [-1,1] range
    if (pol_df['p_polarity'].abs()>1).any():
        tmp_ind=pol_df['p_polarity'].abs().argmax()
        raise ValueError('P polarity weights should range from -1 to 1. Example problem polarity:\n{}'.format(pol_df.loc[[tmp_ind],:]))

    # Discards P polarities with weights below a selected value
    if p_dict['min_polarity_weight']>0:
        drop_flag=pol_df['p_polarity'].abs()<p_dict['min_polarity_weight']
        if any(drop_flag):
            print('Discarded {} polarity measurements with weights less than {}'.format(sum(drop_flag),p_dict['min_polarity_weight']))
            pol_df=pol_df.drop(pol_df[drop_flag].index).reset_index(drop=True)

    # If takeoff and azimuth columns exist, then we'll use precomputed info
    if {'takeoff','azimuth'}.issubset(pol_df.columns):
        if not('takeoff_uncertainty' in pol_df.columns):
            print('*WARNING: No takeoff_uncertainty provided. Assumming 0 deg takeoff uncertainty')
            pol_df['takeoff_uncertainty']=0
        if not('azimuth_uncertainty' in pol_df.columns):
            print('*WARNING: No azimuth_uncertainty provided. Assumming 0 source-receiver azimuthal uncertainty')
            pol_df['azimuth_uncertainty']=0

    # Determines if takeoff/azimuths need (or can) be calculated
    null_index=pol_df[pol_df.loc[:,['takeoff','azimuth','takeoff_uncertainty','azimuth_uncertainty']].isnull().any(axis=1)].index
    if (not(p_dict['vmodel_paths'] and p_dict['stfile'])) | (null_index.empty):
        p_dict['compute_takeoff_azimuth']=False
        if not(null_index.empty):
            print('*WARNING: {} of of the {} polarities have missing takeoff/azimuths. These polarities will be ignored.'.format(len(null_index),len(pol_df)))
            pol_df=pol_df.drop(null_index).reset_index(drop=True)

    return pol_df


def check_model_eq(cat_df,pol_df,p_dict):
    '''
    Checks that earthquake locations will be contained in the model
    '''
    # Ensures the depth range is on the interval p_dict['look_dep'][2]
    if ((p_dict['look_dep'][1]-p_dict['look_dep'][0]) % p_dict['look_dep'][2])!=0:
        new_dep2=p_dict['look_dep'][0]+(np.ceil((p_dict['look_dep'][1]-p_dict['look_dep'][0])/p_dict['look_dep'][2]))*p_dict['look_dep'][2]
        print('Depth range ({}-{} km) for the lookup table is not on the interval {} km. Changing look_dep[1] to {}'.format(p_dict['look_dep'][0],p_dict['look_dep'][1],p_dict['look_dep'][2],new_dep2))
        p_dict['look_dep'][1]=new_dep2

    # Ensures the distance range is on the interval p_dict['look_del'][2]
    if ((p_dict['look_del'][1]-p_dict['look_del'][0]) % p_dict['look_del'][2])!=0:
        new_del2=p_dict['look_del'][0]+(np.ceil((p_dict['look_del'][1]-p_dict['look_del'][0])/p_dict['look_del'][2]))*p_dict['look_del'][2]
        print('Distance range ({}-{} km) for the lookup table is not on the interval {} km. Changing look_del[1] to {}'.
            format(p_dict['look_del'][0],p_dict['look_del'][1],p_dict['look_del'][2],new_del2) )
        p_dict['look_del'][1]=new_del2

    num_source_depth_bins=int((p_dict['look_dep'][1]-p_dict['look_dep'][0])/3+1)
    if num_source_depth_bins>p_dict['nd0']:
        raise ValueError('Given the lookup depth range of {}-{}km with interval {}km look_dep), the {} source depth bins needed exceeds the maximum number of source depth bins {} (nd0).'.format(p_dict['look_dep'][0],p_dict['look_dep'][1],p_dict['look_dep'][2],num_source_depth_bins,p_dict['nd0']))

    if p_dict['delmin']<=0:
        p_dict['delmin']=p_dict['look_del'][0]
    # If the maximum source-receiver distance is not specified, sets it equal to the velocity model extent
    if p_dict['delmax']<=0:
        p_dict['delmax']=p_dict['look_del'][1]

    # Ensures that the catalog earthquake depths are all contained within the modeled interval.
    if len(cat_df)>0:
        if cat_df['origin_depth_km'].max()>p_dict['look_dep'][1]:
            prob_ind=cat_df['origin_depth_km'].argmax()
            if p_dict['allow_hypocenters_outside_table']:
                print('Earthquake depth ({} km) for event {} is greater than maximum lookup table depth ({} km). Setting hypocentral depth to {} km.'.
                    format(cat_df['origin_depth_km'].max(),cat_df.loc[prob_ind,'event_id'],p_dict['look_dep'][1],p_dict['look_dep'][1]))
                cat_df.loc[cat_df['origin_depth_km']>p_dict['look_dep'][1],'origin_depth_km']=p_dict['look_dep'][1]
            else:
                raise ValueError('Earthquake depth ({} km) for event {} is greater than maximum lookup table depth ({} km).'.
                    format(cat_df['origin_depth_km'].max(),cat_df.loc[prob_ind,'event_id'],p_dict['look_dep'][1]))
        if cat_df['origin_depth_km'].min()<p_dict['look_dep'][0]:
            prob_ind=cat_df['origin_depth_km'].argmin()
            if p_dict['allow_hypocenters_outside_table']:
                print('WARNING: Earthquake depth ({} km) for event {} is less than minimum lookup table depth ({} km). Setting hypocentral depth to {} km.'.
                    format(cat_df['origin_depth_km'].min(),cat_df.loc[prob_ind,'event_id'],p_dict['look_dep'][0],p_dict['look_dep'][0]))
                cat_df.loc[cat_df['origin_depth_km']<p_dict['look_dep'][0],'origin_depth_km']=p_dict['look_dep'][0]
            else:
                raise ValueError('Earthquake depth ({} km) for event {} is less than maximum lookup table depth ({} km).'.
                    format(cat_df['origin_depth_km'].min(),cat_df.loc[prob_ind,'event_id'],p_dict['look_dep'][0]))
    if 'origin_depth_km' in pol_df.columns:
        if pol_df['origin_depth_km'].max()>p_dict['look_dep'][1]:
            if p_dict['allow_hypocenters_outside_table']:
                print('WARNING: Earthquake depth ({}) in polarity record is greater than the maximum lookup table depth ({} km). Setting hypocentral depth to {} km.'.
                    format(pol_df['origin_depth_km'].max(),p_dict['look_dep'][1],p_dict['look_dep'][1]))
                pol_df.loc[pol_df['origin_depth_km']>p_dict['look_dep'][1],'origin_depth_km']=p_dict['look_dep'][1]
            else:
                raise ValueError('Earthquake depth ({}) in polarity record is greater than the maximum lookup table depth ({} km)'.
                    format(pol_df['origin_depth_km'].max(),p_dict['look_dep'][1]))
        if pol_df['origin_depth_km'].min()<p_dict['look_dep'][0]:
            if p_dict['allow_hypocenters_outside_table']:
                print('WARNING: Earthquake depth ({}) in polarity record is less than the minimum lookup table depth ({} km). Setting hypocentral depth to {} km.'.
                    format(pol_df['origin_depth_km'].min(),p_dict['look_dep'][0],p_dict['look_dep'][0]))
                pol_df.loc[pol_df['origin_depth_km']<p_dict['look_dep'][0],'origin_depth_km']=p_dict['look_dep'][0]
            else:
                raise ValueError('Earthquake depth ({}) in polarity record is less than the maximum lookup table depth ({} km)'.
                    format(pol_df['origin_depth_km'].min(),p_dict['look_dep'][0]))

    # Warns the user if a permiated earthquake depth could fall outside of the modeled interval.
    if len(cat_df)>0:
        max_eq_perm_depth=(cat_df['origin_depth_km']+cat_df['vert_uncert_km']).max()
        if max_eq_perm_depth>p_dict['look_dep'][1]:
            prob_ind=np.argmax(max_eq_perm_depth)
            print(('*WARNING: Given the earthquake depths and vertical uncertainties, permiated hypocenters could fall '+
            'outside of the maximum lookup table depth of {} km. Example event: max permiated depth of {} km for event {}.'+
            'You may wish to increase the depth of the lookup table.').format(p_dict['look_dep'][1],max_eq_perm_depth,cat_df.loc[prob_ind,'event_id']))

    return pol_df,p_dict
