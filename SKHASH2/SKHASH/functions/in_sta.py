'''
Functions for reading station files and station reversals.
'''

# Standard libraries
import time

# External libraries
import numpy as np
import pandas as pd


def read_station_file(pol_df,p_dict):
    '''
    Reads input station file.
    '''
    if (p_dict['input_format']=='skhash'):
        station_df=read_skhash_station_file(p_dict['stfile'],p_dict['merge_on'])
    elif (p_dict['input_format']=='hash2') | (p_dict['input_format']=='hash3'):
        station_df=read_hash3_station_file(p_dict['stfile'])
    elif (p_dict['input_format']=='hash4'):
        station_df=read_hash4_station_file(p_dict['stfile'])
    else:
        raise ValueError('Unknown station metadata file format ({}).'.format(p_dict['input_format']))

    station_df[p_dict['merge_on']]=station_df[p_dict['merge_on']].astype(str)
    station_df['sta_code']=station_df[p_dict['merge_on']].agg('.'.join, axis=1)

    # Ensures the required columns are present
    if not({'sta_code','station_lat','station_lon','station_depth_km'}).issubset(station_df.columns):
        raise ValueError('The station dataframe does not contain the required columns:\n\t[sta_code, station_lat, station_lon, station_depth_km]')

    if p_dict['require_temporal_match']:
        if not({'start_time','end_time'}).issubset(station_df.columns):
            raise ValueError('When require_temporal_match==True, the station dataframe must contain the additional columns:\n\t[start_time, end_time]')
        # Converts the start/end time strings into timestamps.
        station_df['start_time']=pd.to_datetime(station_df['start_time'],infer_datetime_format=True,errors='coerce')
        station_df['start_time']=station_df['start_time'].fillna(pd.Timestamp('1700-01-01'))
        station_df['end_time']=pd.to_datetime(station_df['end_time'],infer_datetime_format=True,errors='coerce')
        station_df['end_time']=station_df['end_time'].fillna(pd.Timestamp('2200-01-01'))

        # Selects only the desired columns
        station_df=station_df.filter(['sta_code','station_lat','station_lon','station_depth_km','start_time','end_time'])
    else:
        # Selects only the desired columns
        station_df=station_df.filter(['sta_code','station_lat','station_lon','station_depth_km'])

    # Selects only the station metadata that has polarity or S/P amp information
    station_df=station_df.loc[station_df['sta_code'].isin(pol_df['sta_code'].unique()),:].reset_index(drop=True)

    # Drops duplicate records
    station_df=station_df.drop_duplicates().reset_index(drop=True)

    return station_df


def read_skhash_station_file(stfile,merge_on):
    '''
    Reads input station file following the SKHASH format.
    '''
    station_df=pd.read_csv(stfile)
    required_column_names=['latitude','longitude','elevation']
    if not(set(required_column_names).issubset(station_df.columns)):
        raise ValueError('When using the SKHASH input format for station metadata files, the following header names MUST be provided:\n{}'.format(required_column_names))
    for req_col in merge_on:
        if not(req_col in station_df.columns):
            raise ValueError('When using SKHASH input format and require_{}_match==True, the station file (stfile: {}) must contain the column "{}".'.format(req_col,stfile,req_col))

    station_df=station_df.rename(columns={'latitude':'station_lat','longitude':'station_lon','elevation':'station_elevation_m'})
    # station_df=pd.read_csv(stfile, header=None,names=['station','network','location','channel','station_lat','station_lon','station_elevation_m','start_time','end_time'], delimiter=r"\s+")
    station_df['station_depth_km']=station_df['station_elevation_m']/1000*-1 # Converts station elevation (in m) to station depth (in km)
    station_df=station_df.drop(columns='station_elevation_m')

    return station_df


def read_hash3_station_file(stfile):
    '''
    Reads input station file following the (old) SCEDC format used for HASH Driver 2, 3, and 5.
    '''
    station_df=pd.read_fwf(stfile,colspecs=[[0,4],[5,8],[41,50],[51,61],[62,67],[68,78],[79,89],[90,92]],names=['station','channel','station_lat','station_lon','station_elevation_m','start_time','end_time','network'])
    station_df['network']=station_df['network'].str.strip()
    station_df['station']=station_df['station'].str.strip()
    station_df['channel']=station_df['channel'].str.strip()
    station_df['station_depth_km']=station_df['station_elevation_m']/1000*-1
    station_df.drop(columns='station_elevation_m',inplace=True)
    station_df=station_df.loc[:,['station','network','channel','station_lat','station_lon','station_depth_km','start_time','end_time']]

    return station_df


def read_hash4_station_file(stfile):
    '''
    Reads input station file following the (newer) SCEDC format used for HASH Driver 4.
    '''
    station_df=pd.read_fwf(stfile,colspecs=[[0,2],[4,9],[10,13],[16,18],[51,60],[60,71],[71,77],[77,88],[88,99]],names=['network','station','channel','location','station_lat','station_lon','station_elevation_m','start_time','end_time'])
    station_df['network']=station_df['network'].str.strip()
    station_df['station']=station_df['station'].str.strip()
    station_df['channel']=station_df['channel'].str.strip()
    station_df['station_depth_km']=station_df['station_elevation_m']/1000*-1
    station_df.drop(columns='station_elevation_m',inplace=True)
    station_df=station_df.loc[:,['station','network','channel','station_lat','station_lon','station_depth_km','start_time','end_time']]

    return station_df


def apply_station_locations(pol_df,station_df,p_dict):
    '''
    Mergest the station locations to the polarity information.
    '''

    # Checks to see if there is missing metadata for the polarities
    missing_metadata_index=pol_df[~(pol_df['sta_code'].isin(station_df['sta_code']))].index
    if len(missing_metadata_index)>0:
        missing_sta_df=pol_df.loc[missing_metadata_index,'sta_code'].drop_duplicates().reset_index(drop=True)
        num_missing_sta=len(missing_sta_df)
        if p_dict['ignore_missing_metadata']:
            print('WARNING: Missing metadata for {} measurements. These measurements will be ignored. Example missing stations:\n{}'.format(num_missing_sta,missing_sta_df.head(10).to_string()))
            pol_df=pol_df.drop(missing_metadata_index).reset_index(drop=True)
        else:
            raise ValueError('Missing metadata for {} measurements. Example missing stations:\n{}'.format(num_missing_sta,missing_sta_df.head(10).to_string()))

    if p_dict['require_temporal_match']:
        time_merge_event_sta=time.time()

        pol_df=pol_df.reset_index(drop=True)
        station_df=station_df.reset_index(drop=True)

        # Merges measurement records and station metadata considering receiver start/end times.
        # This is slow and memory intensive, and it could be greatly improved!
        i, j = np.where((pol_df.origin_DateTime.values[:, None] >= station_df.start_time.values) &
                        (pol_df.origin_DateTime.values[:, None] <= station_df.end_time.values) &
                        (pol_df.sta_code.values[:,None]==station_df.sta_code.values))

        found_pol_df=pd.concat([pol_df.loc[i, :].reset_index(drop=True),
                            station_df.loc[j, ['station_lat','station_lon','station_depth_km','start_time','end_time']].reset_index(drop=True)], axis=1)
        missing_pol_df=pol_df.iloc[~np.in1d(np.arange(len(pol_df)), np.unique(i)),:]

        if len(missing_pol_df)>0:
            print(('{} measurements did not match station metadata considering station start/end times.\n'+
                    'Ignoring start/end times for these will be ignored, and the first matching metadata record will be used instead.'+
                    'Example record:\n{}\n').format(len(missing_pol_df),missing_pol_df.loc[:,['event_id','sta_code','origin_DateTime']].head(10)))

        # If there are overlapping metadata windows, more than one record could be applied. This selects the first matching station.
        found_pol_df=found_pol_df.drop_duplicates(subset=['sta_code','event_id','source']).reset_index(drop=True)

        # Makes sure the number records with found and missing metadata match the original number of records
        if (len(found_pol_df)+len(missing_pol_df))!=len(pol_df):
            raise ValueError('Unexpected error while matching station start/times with polarities. The number of found+missing polarities do not match the original number of polarities.')

        # Ignores the start/end times and appends the first matching metadata record for each measurement
        missing_pol_df=missing_pol_df.merge(station_df.drop_duplicates(subset=['sta_code']),on='sta_code',how='left')

        # Concats the found/previously missing records back together
        pol_df=pd.concat([found_pol_df,missing_pol_df]).sort_values(by=['event_id','sta_code']).reset_index(drop=True)

        print('Measurements merged with station metadata considering start/end times. Runtime: {:.2f} sec'.format(time.time()-time_merge_event_sta))

    else:

        # Looks to see if there are multiple locations for the same metadata selection. If so, the first listed location is used
        duplicate_station_index=station_df[station_df.duplicated(subset='sta_code')].index
        if not(duplicate_station_index.empty):
            print(('*WARNING: {} station records have duplicate locations for your selected metadata attributes: {}\n'+
                            '\tThe first station location will be used. Example problematic station info:\n{}\n').
                            format(len(duplicate_station_index),p_dict['merge_on'],station_df.loc[duplicate_station_index,:].reset_index(drop=True).head(10).to_string()))
            station_df=station_df.drop(duplicate_station_index).reset_index(drop=True)

        pol_df=pol_df.merge(station_df[['sta_code','station_lat','station_lon','station_depth_km']],on='sta_code',how='left')

    missing_metadata_index=pol_df[pd.isnull(pol_df.loc[:,['station_lat','station_lon','station_depth_km']]).any(axis=1)].index
    if len(missing_metadata_index)>0:
        missing_sta_df=pol_df.loc[missing_metadata_index,'sta_code'].drop_duplicates().reset_index(drop=True)
        raise ValueError('Error: Missing metadata for {} stations:\n{}'.format(len(missing_sta_df),missing_sta_df.to_string()))

    # Calculates source-receiver distances
    aspect=np.cos(np.deg2rad(pol_df['origin_lat'].values))
    dx=(pol_df['station_lon'].values-pol_df['origin_lon'].values)*111.2*aspect
    dy=(pol_df['station_lat'].values-pol_df['origin_lat'].values)*111.2
    pol_df['sr_dist_km']=np.sqrt(dx**2+dy**2)

    return pol_df


def read_reverse_file(p_dict):
    '''
    Reads input station reversal file.
    '''
    if (p_dict['input_format']=='skhash'):
        pol_reverse_df=read_reverse_skhash_file(p_dict['plfile'],p_dict['merge_on'])
    elif (p_dict['input_format'][:-1]=='hash'):
        pol_reverse_df=read_reverse_hash_file(p_dict['plfile'])
    else:
        raise ValueError('Unknown station polarity reversal file format ({}).'.format(p_dict['input_format']))
    return pol_reverse_df


def read_reverse_skhash_file(plfile,merge_on):
    '''
    Reads input station reversal file using the SKHASH format.
    '''
    try:
        pol_reverse_df=pd.read_csv(plfile,skipinitialspace=True,comment='#')
    except pd.errors.EmptyDataError:
        print('The polarity reversal file (plfile: {}) is empty.'.format(plfile))
        pol_reverse_df=pd.DataFrame(columns=np.hstack([merge_on,'start_time','end_time']))

    required_column_names=['start_time','end_time']
    if not(set(required_column_names).issubset(pol_reverse_df.columns)):
        raise ValueError('When using the skhash input format for station metadata files, the following header names MUST be provided:\n{}'.format(required_column_names))
    for req_col in merge_on:
        if not(req_col in pol_reverse_df.columns):
            raise ValueError('When using SKHASH input format and require_{}_match==True, the station polarity reversal file (plfile: {}) must contain the column "require_{}_match".'.format(req_col,plfile,req_col))

    pol_reverse_df.loc[(pol_reverse_df['start_time'].astype(str)=='0'),'start_time']='1900-01-01'
    pol_reverse_df.loc[(pol_reverse_df['end_time'].astype(str)=='0'),'end_time']='2200-01-01'

    pol_reverse_df['start_time']=pd.to_datetime(pol_reverse_df['start_time'])
    pol_reverse_df['end_time']=pd.to_datetime(pol_reverse_df['end_time'])
    return pol_reverse_df


def read_reverse_hash_file(plfile):
    '''
    Reads input station reversal file using the HASH format.
    '''
    pol_reverse_df=pd.read_csv(plfile,sep='\s+',names=['station','start_time_int','end_time_int'])
    pol_reverse_df.loc[pol_reverse_df['start_time_int']==0,'start_time_int']=19000101
    pol_reverse_df.loc[pol_reverse_df['end_time_int']==0,'end_time_int']=22000101
    pol_reverse_df['start_time']=pd.to_datetime(pol_reverse_df['start_time_int'],format='%Y%m%d')
    pol_reverse_df['end_time']=pd.to_datetime(pol_reverse_df['end_time_int'],format='%Y%m%d')
    pol_reverse_df=pol_reverse_df.drop(columns=['start_time_int','end_time_int'])
    return pol_reverse_df


def reverse_polarities(pol_df,pol_reverse_df,p_dict):
    '''
    Flips the polarities for reversed stations.
    '''
    if (not(pol_df.empty)) & (not(pol_reverse_df.empty)):
        if (p_dict['input_format']=='skhash'):

            if not('origin_DateTime' in pol_df.columns):
                print('Applying station reversals ignoring the origin time.')

            # Computes "sta codes" by combining the desired metadata (network, station, location, channel) codes
            pol_reverse_df['sta_code']=pol_reverse_df[p_dict['merge_on']].agg('.'.join, axis=1)

            # Only considers reversed stations that have a polarity measurement
            pol_reverse_df=pol_reverse_df.loc[pol_reverse_df['sta_code'].isin(pol_df['sta_code']),:].reset_index(drop=True)

            # Determines the indecies to be reversed
            potential_reverse_ind=pol_df.index[pol_df['sta_code'].isin(pol_reverse_df['sta_code'])]
            potential_reverse_sta_code=pol_df.loc[potential_reverse_ind,'sta_code'].unique()
            group_potential_reverse_pol_df=pol_df.loc[potential_reverse_ind,:].groupby('sta_code')
            reverse_indicies=[]
            for sta_code in potential_reverse_sta_code:
                tmp_df=group_potential_reverse_pol_df.get_group(sta_code)
                sta_pol_reverse_df=pol_reverse_df.loc[pol_reverse_df['sta_code']==sta_code,:].reset_index(drop=True)

                if 'origin_DateTime' in tmp_df.columns:
                    for tmp_x in range(len(sta_pol_reverse_df)):
                        reverse_indicies.append(tmp_df.index[ ((tmp_df['origin_DateTime']>sta_pol_reverse_df.loc[tmp_x,'start_time']) & (tmp_df['origin_DateTime']<sta_pol_reverse_df.loc[tmp_x,'end_time'])) ])
                else:
                    reverse_indicies.append(tmp_df.index)
            if len(reverse_indicies)>0:
                reverse_indicies=np.unique(np.hstack(reverse_indicies))

            pol_df.loc[reverse_indicies,'p_polarity']=pol_df.loc[reverse_indicies,'p_polarity']*-1

        elif (p_dict['input_format'][:-1]=='hash'):
            # Only considers reversed stations that have a polarity measurement
            pol_reverse_df=pol_reverse_df.loc[pol_reverse_df['station'].isin(pol_df['station']),:].reset_index(drop=True)

            # Determines the indecies to be reversed
            potential_reverse_ind=pol_df.index[pol_df['station'].isin(pol_reverse_df['station'])]
            potential_reverse_sta=pol_df.loc[potential_reverse_ind,'station'].unique()
            group_potential_reverse_pol_df=pol_df.loc[potential_reverse_ind,:].groupby('station')
            reverse_indicies=[]
            for sta in potential_reverse_sta:
                tmp_df=group_potential_reverse_pol_df.get_group(sta)
                sta_pol_reverse_df=pol_reverse_df.loc[pol_reverse_df['station']==sta,:].reset_index(drop=True)
                for tmp_x in range(len(sta_pol_reverse_df)):
                    reverse_indicies.append(tmp_df.index[ ((tmp_df['origin_DateTime']>sta_pol_reverse_df.loc[tmp_x,'start_time']) & (tmp_df['origin_DateTime']<sta_pol_reverse_df.loc[tmp_x,'end_time'])) ])
            reverse_indicies=np.unique(np.hstack(reverse_indicies))

            pol_df.loc[reverse_indicies,'p_polarity']=pol_df.loc[reverse_indicies,'p_polarity']*-1

        else:
            raise ValueError('Unknown station polarity reversal file format ({}).'.format(p_dict['input_format']))

    return pol_df
