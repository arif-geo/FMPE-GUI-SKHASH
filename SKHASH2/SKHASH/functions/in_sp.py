'''
Functions for reading and computing S/P amplitude ratios.
'''

# External libraries
import numpy as np
import pandas as pd


def read_amp_corr_files(ampfile_in,p_dict):
    '''
    Reads the S/P ratio file (ampfile) and applies the corrections (corfile)
    '''
    # Reads S/P file and calculates ratios
    spamp_df=read_amp_file(ampfile_in,p_dict['input_format'],p_dict['merge_on'],p_dict['allow_duplicate_stations'],p_dict['ratmin'],p_dict['min_sp'],p_dict['max_sp'])

    if spamp_df.empty:
        p_dict['ampfile']=''
    else:
        # Reads station corrections and applies them to the S/P measurements
        if p_dict['corfile']:
            stacor_df=read_sta_corr(p_dict)
            spamp_df=apply_sta_correction(spamp_df,stacor_df)

    return spamp_df,p_dict


def read_sta_corr(p_dict):
    '''
    Reads file of station corrections. Returns a data frame containing the
    "sta codes" and their corresponding correction.
    '''
    if p_dict['input_format']=='skhash':
        stacor_df=read_skhash_sta_corr(p_dict['corfile'],p_dict['merge_on'])
    elif p_dict['input_format']=='hash3':
        stacor_df=read_hash3_sta_corr(p_dict['corfile'])
    else:
        raise ValueError('Unknown S/P amplitude station correction file file format ({}).'.format(p_dict['input_format']))

    # Creates "sta codes" for the corrections
    stacor_df['sta_code']=stacor_df[p_dict['merge_on']].agg('.'.join, axis=1)

    # Checks for any duplicate station correction entries
    stacor_dup_flag=stacor_df['sta_code'].duplicated()
    if stacor_dup_flag.any():
        print('*WARNING: There are {} duplicate station corrections. Only the first occurrence of each will be kept. Example duplicate:\n{}'.
            format(np.sum(stacor_dup_flag),stacor_df.loc[[np.argmax(stacor_dup_flag)],:].to_string()))
        stacor_df=stacor_df.drop(stacor_df[stacor_dup_flag].index).reset_index(drop=True)

    stacor_df=stacor_df[['sta_code','sta_correction']]

    return stacor_df


def read_hash3_sta_corr(corfile):
    '''
    Reads file of station corrections using the HASH3 format
    '''
    stacor_df=pd.read_csv(corfile,sep='\s+',names=['station','channel','network','sta_correction'])
    return stacor_df


def read_skhash_sta_corr(corfile,merge_on):
    '''
    Reads file of station corrections using the SKHASH format
    '''
    consider_cols=['network','station','location','channel','sta_correction']
    try:
        stacor_df=pd.read_csv(corfile,skipinitialspace=True,comment='#',usecols=lambda x: x in consider_cols)
    except pd.errors.EmptyDataError:
        print('*Warning: corfile ({}) is empty.'.format(corfile))
        stacor_df=[]
    if not('sta_correction' in stacor_df.columns):
        raise ValueError('When using SKHASH input format, the corfile ({}) must contain the column "sta_correction".'.format(corfile))
    for req_col in merge_on:
        if not(req_col in stacor_df.columns):
            raise ValueError('When using SKHASH input format and require_{}_match==True, the corfile ({}) must contain the column "{}".'.format(req_col,corfile,req_col))
    return stacor_df


def read_amp_file(ampfile,input_format,merge_on,allow_duplicate_stations,ratmin,min_sp,max_sp):
    '''
    Reads file of S/P ratios. Returns a data frame containing the
    "sta codes" and their corresponding S/P ratios.
    '''
    if input_format=='skhash':
        spamp_df=read_skhash_amp_file(ampfile,merge_on)
    elif input_format=='hash3':
        spamp_df=read_hash3_amp_file(ampfile)
    else:
        raise ValueError('Unknown S/P amplitude file file format ({}).'.format(input_format))

    spamp_df['event_id']=spamp_df['event_id'].astype(str)
    if 'event_id2' in spamp_df.columns:
        spamp_df['event_id2']=spamp_df['event_id2'].astype(str)

    if ('sp_ratio' in spamp_df.columns):
        spamp_df['sp_ratio']=spamp_df['sp_ratio'].abs()
    else:
        # spamp_df.loc[:,['noise_p','noise_s','amp_p','amp_s']]=spamp_df.loc[:,['noise_p','noise_s','amp_p','amp_s']].abs()
        spamp_df[['noise_p','noise_s','amp_p','amp_s']]=spamp_df[['noise_p','noise_s','amp_p','amp_s']].abs()

        spamp_df=spamp_df.loc[((spamp_df['amp_p']/spamp_df['noise_p'])>=ratmin) &
                                ((spamp_df['amp_s']/spamp_df['noise_s'])>=ratmin),:].reset_index(drop=True)
        spamp_df['sp_ratio']=(spamp_df['amp_s']/spamp_df['amp_p']).round(2)
        spamp_df=spamp_df.drop(columns=['noise_p','noise_s','amp_p','amp_s'])

    spamp_df=spamp_df.drop(spamp_df[spamp_df['sp_ratio']==0].index).reset_index(drop=True)

    # Sets the minimum and maximum S/P ratios
    if min_sp>0:
        spamp_df.loc[spamp_df.sp_ratio<min_sp,'sp_ratio']=min_sp
    if max_sp>0:
        spamp_df.loc[spamp_df.sp_ratio>max_sp,'sp_ratio']=max_sp

    spamp_df['sp_ratio']=np.log10(spamp_df['sp_ratio'])
    spamp_df['sta_code']=spamp_df[merge_on].agg('.'.join, axis=1)

    if (not(allow_duplicate_stations)) and (len(spamp_df)>0):
        # Checks for any duplicate S/P entries
        spamp_dup_flag=spamp_df[['sta_code','event_id']].duplicated()
        if spamp_dup_flag.any():
            print('*WARNING: There are {} duplicate S/P amplitudes. Only the first occurrence of each will be kept. Example duplicate:\n{}'.
                format(np.sum(spamp_dup_flag),spamp_df.loc[[np.argmax(spamp_dup_flag)],:].to_string()))
            spamp_df=spamp_df.drop(spamp_df[spamp_dup_flag].index).reset_index(drop=True)

    spamp_df=spamp_df.filter(['event_id','event_id2','sta_code','sp_ratio','origin_lat','origin_lon','origin_depth_km','takeoff','takeoff_uncertainty','azimuth','azimuth_uncertainty'])

    return spamp_df


def read_skhash_amp_file(ampfile,merge_on):
    '''
    Reads file of S/P ratios using the SKHASH format
    '''
    consider_cols=['event_id','event_id2','network','station','location','channel','noise_p','noise_s','amp_p','amp_s','sp_ratio','origin_latitude','origin_longitude','origin_depth_km','takeoff','takeoff_uncertainty','azimuth','azimuth_uncertainty']
    try:
        spamp_df=pd.read_csv(ampfile,skipinitialspace=True,usecols=lambda x: x in consider_cols)
    except pd.errors.EmptyDataError:
        print('*Warning: ampfile ({}) is empty.'.format(ampfile))
        spamp_df=[]
    if not('event_id' in spamp_df.columns):
        raise ValueError('When using SKHASH input format, the ampfile ({}) must contain the column "event_id".'.format(ampfile))
    if (not({'noise_p','noise_s','amp_p','amp_s'}.issubset(spamp_df.columns))) and (not('sp_ratio' in spamp_df.columns)):
        raise ValueError('When using SKHASH input format, the ampfile ({}) must either contain one of the following sets of column names:\n\t1. noise_p,noise_s,amp_p,amp_s\n\t2. sp_ratio'.format(ampfile))

    for req_col in merge_on:
        if not(req_col in spamp_df.columns):
            raise ValueError('When using SKHASH input format and require_{}_match==True, the ampfile ({}) must contain the column "{}".'.format(req_col,ampfile,req_col))

    spamp_df=spamp_df.rename(columns={'origin_latitude':'origin_lat','origin_longitude':'origin_lon'})

    return spamp_df


def read_hash3_amp_file(ampfile):
    '''
    Reads file of S/P amplitude ratios using the HASH3 format
    '''
    file1 = open(ampfile, 'r')
    event_amplitude_list=[]
    sta=[];net=[];loc=[];cha=[];noise_p=[];noise_s=[];amp_p=[];amp_s=[]
    event_id=''
    for line in file1:
        split_line=line.split()
        if len(split_line)==2: # Event id line
            if len(sta)>0:
                tmp_df=pd.DataFrame(data=np.asarray(np.vstack((sta,net,loc,cha,noise_p,noise_s,amp_p,amp_s))).T,
                                    columns=['station','network','location','channel','noise_p','noise_s','amp_p','amp_s'])
                tmp_df['event_id']=event_id
                tmp_df=tmp_df.astype({'noise_p':float,'noise_s':float,'amp_p':float,'amp_s':float})
                event_amplitude_list.append(tmp_df)
                sta=[];net=[];loc=[];cha=[];noise_p=[];noise_s=[];amp_p=[];amp_s=[]
            event_id=split_line[0]
        else: # Amplitude line
            sta.append(split_line[0])
            cha.append(split_line[1])
            net.append(split_line[2])
            noise_p.append(split_line[5])
            noise_s.append(split_line[6])
            amp_p.append(split_line[7])
            amp_s.append(split_line[8])
            loc.append('--')
    if len(sta)>0:
        tmp_df=pd.DataFrame(data=np.asarray(np.vstack((sta,net,loc,cha,noise_p,noise_s,amp_p,amp_s))).T,
                            columns=['station','network','location','channel','noise_p','noise_s','amp_p','amp_s'])
        tmp_df['event_id']=event_id
        tmp_df=tmp_df.astype({'noise_p':float,'noise_s':float,'amp_p':float,'amp_s':float})
        event_amplitude_list.append(tmp_df)
    file1.close()

    if len(event_amplitude_list)>0:
        spamp_df=pd.concat(event_amplitude_list).reset_index(drop=True)
    else:
        spamp_df=pd.DataFrame(columns=['station','network','location','channel','noise_p','noise_s','amp_p','amp_s'])

    return spamp_df


def apply_sta_correction(spamp_df,stacor_df):
    '''
    Applies station corrections to S/P ratios.
    '''
    spamp_df=spamp_df.merge(stacor_df,on='sta_code',how='left')
    spamp_df.loc[spamp_df[pd.isnull(spamp_df['sta_correction'])].index,'sta_correction']=0
    spamp_df['sp_ratio']-=spamp_df['sta_correction']
    spamp_df=spamp_df.drop(columns='sta_correction')

    return spamp_df
