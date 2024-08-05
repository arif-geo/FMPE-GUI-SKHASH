'''
Functions for reading P-polarity inputs.
'''

# External libraries
import numpy as np
import pandas as pd


def read_polarity_file(fpfile,input_format,merge_on):
    if (input_format=='skhash'):
        tmp_pol_df=read_skhash_polarity_file(fpfile,merge_on)
        cat_df=[]
    elif (input_format=='ncsn') or (input_format=='hypoinverse'):
        cat_df,tmp_pol_df=read_ncsn_polarity_file(fpfile)
    elif input_format=='hash1':
        cat_df,tmp_pol_df=read_hash1_polarity_file(fpfile)
    elif (input_format=='hash2') | (input_format=='hash3') | (input_format=='hash5'):
        cat_df,tmp_pol_df=read_hash3_polarity_file(fpfile)
    elif (input_format=='hash4'):
        cat_df,tmp_pol_df=read_hash4_polarity_file(fpfile)
    else:
        raise ValueError('Unknown p-wave first motion polarity file format (fpfile: {}).'.format(input_format))

    # Ensures column data formats are the expected type
    tmp_pol_df['event_id']=tmp_pol_df['event_id'].astype(str)
    tmp_pol_df['p_polarity']=tmp_pol_df['p_polarity'].astype(float)
    for merge_on_it in merge_on:
        tmp_pol_df[merge_on_it]=tmp_pol_df[merge_on_it].astype(str)
    if 'event_id2' in tmp_pol_df.columns:
        tmp_pol_df['event_id2']=tmp_pol_df['event_id2'].astype(str)
    if 'sr_dist_km' in tmp_pol_df.columns:
        tmp_pol_df['sr_dist_km']=tmp_pol_df['sr_dist_km'].astype(float)
    if 'takeoff' in tmp_pol_df.columns:
        tmp_pol_df['takeoff']=tmp_pol_df['takeoff'].astype(float)
    if 'azimuth' in tmp_pol_df.columns:
        tmp_pol_df['azimuth']=tmp_pol_df['azimuth'].astype(float)
    if 'takeoff_uncertainty' in tmp_pol_df.columns:
        tmp_pol_df['takeoff_uncertainty']=tmp_pol_df['takeoff_uncertainty'].astype(float)
    if 'azimuth_uncertainty' in tmp_pol_df.columns:
        tmp_pol_df['azimuth_uncertainty']=tmp_pol_df['azimuth_uncertainty'].astype(float)

    # Computes "sta codes" by combining the desired metadata (network, station, location, channel) codes
    tmp_pol_df['sta_code']=tmp_pol_df[merge_on].agg('.'.join, axis=1)

    # List of desired columns. We have to keep station codes for HASH simulps and polarity reversal compatability
    require_col=['event_id','sta_code','p_polarity','station','sr_dist_km','takeoff','azimuth','takeoff_uncertainty','azimuth_uncertainty']
    additional_col=['event_id2','origin_lat','origin_lon','origin_depth_km']

    # Selects only the desired columns.
    tmp_pol_df=tmp_pol_df.filter(np.hstack((require_col,additional_col)))

    # Ensures all desired columns exist
    for col in require_col:
        if col not in tmp_pol_df.columns:
            tmp_pol_df[col]=np.nan

    return cat_df,tmp_pol_df


def read_skhash_polarity_file(fpfile,merge_on):
    '''
    Reads polarity file using the SKHASH format.
    '''
    consider_cols=['event_id','event_id2','network','station','location','channel','p_polarity','takeoff','takeoff_uncertainty','azimuth','azimuth_uncertainty','sr_dist_km','origin_latitude','origin_longitude','origin_depth_km','horz_uncert_km','vert_uncert_km']
    try:
        tmp_pol_df=pd.read_csv(fpfile,skipinitialspace=True,comment='#',usecols=lambda x: x in consider_cols)
    except pd.errors.EmptyDataError:
        raise ValueError('fpfile ({}) is empty.'.format(fpfile))
    if not({'event_id', 'p_polarity'}.issubset(tmp_pol_df.columns)):
        raise ValueError('When using SKHASH input format, the fpfile ({}) must contain the columns "event_id" and "p_polarity".'.format(fpfile))
    for req_col in merge_on:
        if not(req_col in tmp_pol_df.columns):
            raise ValueError('When using SKHASH input format and require_{}_match==True, the fpfile ({}) must contain the column "{}".'.format(req_col,fpfile,req_col))

    tmp_pol_df=tmp_pol_df.rename(columns={'origin_latitude':'origin_lat','origin_longitude':'origin_lon'})

    return tmp_pol_df


def dm2dd(deg, min, hemisphere):
    '''
    Converts degrees-minutes to decimal degrees
    '''
    dd = float(deg) + float(min)/60
    if hemisphere == 'W' or hemisphere == 'S':
        dd *= -1
    return dd


def read_ncsn_polarity_file(fpfile,use_weight_code=True):
    '''
    Reads a USGS NCSN catalog + phase polarity file in hypoinverse format.
    Polarities can be obtained from the NCEDC (https://ncedc.org/ncedc/catalog-search.html)
    and selecting "NCSN catalog + Phase in Hypoinverse format".

    If use_weight_code==False:
        Gives Impulsive (I) manual picks a weight of p_weight_I (default 1.0)
        Gives Emmergent (E) manual picks a weight of p_weight_E (default 0.5)
        If no I/E flag is provided, assumes E.
    If use_weight_code==True:
        pick_df.loc[pick_df['p_weight_code']==0,'p_polarity']=1
        pick_df.loc[pick_df['p_weight_code']==1,'p_polarity']=0.5
        pick_df.loc[pick_df['p_weight_code']==2,'p_polarity']=0.2
        pick_df.loc[pick_df['p_weight_code']==3,'p_polarity']=0.1
        (see https://ncedc.org/pub/doc/ncsn/shadow2000.pdf) for the weight of the polarity.
    '''
    p_weight_I=1
    p_weight_E=0.5
    col_ind_name=np.asarray([
        [0,4,'year'],
        [4,6,'month'],
        [6,8,'day'],
        [8,10,'hour'],
        [10,12,'minute'],
        [12,16,'second'],
        [16,18,'lat_deg'],
        [18,19,'latNS'],
        [19,23,'lat_min'],
        [23,26,'lon_deg'],
        [26,27,'lonEW'],
        [27,31,'lon_min'],
        [31,36,'depth'],
        [85,89,'x_error_km'],
        [89,93,'z_error_km'],
        [136,146,'event_id'],
        [147,150,'mag_pref']
    ])
    col_ind=col_ind_name[:,0:2].astype(int)
    col_name=col_ind_name[:,2]

    pick_col_ind_name=np.asarray([
        [0,5,'station'],
        [5,7,'network'],
        [9,12,'channel'],
        [13,14,'p_onset'],
        [15,16,'p_first_motion'],
        [16,17,'p_weight_code'],
        [111,113,'location'],
    ])
    pick_col_ind=pick_col_ind_name[:,0:2].astype(int)
    pick_col_name=pick_col_ind_name[:,2]

    header_all=[]
    pick_all=[]
    with open(fpfile, 'r') as file1:
        for line in file1:
            if line[0:3]=='   ': # Footer line
                tmp_pick_df=pd.DataFrame(tmp_pol,columns=pick_col_name)

                # Drops any phase arrivals that do not have a P first motion
                tmp_pick_df=tmp_pick_df.drop(tmp_pick_df.loc[tmp_pick_df['p_first_motion']==' '].index)
                # Adds event_id to the dataframe
                tmp_pick_df['event_id']=str(event_id)

                pick_all.append(tmp_pick_df)
            elif (' ' in line[:10]): # Phase arrival line
                tmp_pol.append([line[i:j] for i,j in zip(pick_col_ind[:,0], pick_col_ind[:,1])])
            else: # Header line
                header_all.append([line[i:j] for i,j in zip(col_ind[:,0], col_ind[:,1])])
                header_all[-1][-2]=header_all[-1][-2].strip()
                event_id=header_all[-1][-2]
                tmp_pol=[]
    cat_df=pd.DataFrame(header_all,columns=col_name)
    cat_df=cat_df.astype({'lat_deg':'int32','lat_min':'int32',
                                'lon_deg':'int32','lon_min':'int32',
                                'depth':'int32','x_error_km':'int32',
                                'z_error_km':'int32','mag_pref':'int32',
                                'year':'int32','month':'int32',
                                'day':'int32','hour':'int32',
                                'minute':'int32','second':'float',
                        })
    cat_df['second']=cat_df['second']/100
    cat_df['origin_DateTime']=pd.to_datetime(cat_df.loc[:,['year','month','day','hour','minute','second']])
    cat_df['origin_lat']=cat_df['lat_deg']+cat_df['lat_min']/6000
    cat_df.loc[cat_df['latNS']!=' ','origin_lat']*=-1
    cat_df['origin_lat']=cat_df['origin_lat'].round(6)
    cat_df['origin_lon']=cat_df['lon_deg']+cat_df['lon_min']/6000
    cat_df.loc[cat_df['lonEW']==' ','origin_lon']*=-1
    cat_df['origin_lon']=cat_df['origin_lon'].round(6)
    cat_df['origin_depth_km']=cat_df['depth']/100
    cat_df['horz_uncert_km']=cat_df['x_error_km']/100
    cat_df['vert_uncert_km']=cat_df['z_error_km']/100
    cat_df['event_mag']=cat_df['mag_pref']/100

    cat_df=cat_df.loc[:,['origin_DateTime', 'origin_lat', 'origin_lon', 'origin_depth_km',
                        'horz_uncert_km', 'vert_uncert_km', 'event_mag', 'event_id']]

    pick_df=pd.concat(pick_all).reset_index(drop=True)
    pick_df['p_polarity']=0

    if use_weight_code:
        pick_df['p_weight_code']=pick_df['p_weight_code'].astype(int)
        pick_df.loc[pick_df['p_weight_code']==0,'p_polarity']=1
        pick_df.loc[pick_df['p_weight_code']==1,'p_polarity']=0.5
        pick_df.loc[pick_df['p_weight_code']==2,'p_polarity']=0.2
        pick_df.loc[pick_df['p_weight_code']==3,'p_polarity']=0.1
    else:
        pick_df.loc[pick_df['p_onset']=='I','p_polarity']=p_weight_I
        pick_df.loc[(pick_df['p_onset']=='E') | (pick_df['p_onset']==' '),'p_polarity']=p_weight_E

    pick_df.loc[pick_df['p_first_motion']=='D','p_polarity']*=-1
    pick_df=pick_df.drop(columns=['p_first_motion','p_weight_code','p_onset'])

    pick_df['station']=pick_df['station'].str.strip()
    pick_df['network']=pick_df['network'].str.strip()
    pick_df['channel']=pick_df['channel'].str.strip()
    pick_df['location']=pick_df['location'].str.strip()

    return cat_df,pick_df


def read_hash1_polarity_file(fpfile,p_weight_I=1.0,p_weight_E=0.5):
    '''
    Reads a polarity file following HASH driver 1 format.
    Gives Impulsive (I) manual picks a weight of p_weight_I (default 1.0)
    Gives Emmergent (E) manual picks a weight of p_weight_E (default 0.5)
    If no I/E flag is provided, assumes E.
    '''
    col_ind_name=np.asarray([
        [0,2,'year'],
        [2,4,'month'],
        [4,6,'day'],
        [6,8,'hour'],
        [8,10,'minute'],
        [10,14,'second'],
        [14,16,'lat_deg'],
        [16,17,'latNS'],
        [17,21,'lat_min'],
        [21,24,'lon_deg'],
        [24,25,'lonEW'],
        [25,29,'lon_min'],
        [29,34,'depth'],
        [34,36,'mag_pref'],
        [80,84,'x_error_km'],
        [84,88,'z_error_km'],
        [122,138,'event_id'],
    ])
    col_ind=col_ind_name[:,0:2].astype(int)
    col_name=col_ind_name[:,2]

    pick_col_ind_name=np.asarray([
        [0,4,'station'],
        [6,7,'p_first_motion'],
        [7,8,'p_weight_code'],
        [58,62,'sr_dist_km'],
        [62,66,'takeoff'],
        [75,78,'azimuth'],
        [79,82,'takeoff_uncertainty'],
        [83,86,'azimuth_uncertainty'],
        [95,98,'channel'],
    ])
    pick_col_ind=pick_col_ind_name[:,0:2].astype(int)
    pick_col_name=pick_col_ind_name[:,2]

    header_all=[]
    pick_all=[]
    with open(fpfile, 'r') as file1:
        is_header_line=True
        for line in file1:
            if is_header_line:
                header_all.append([line[i:j] for i,j in zip(col_ind[:,0], col_ind[:,1])])
                header_all[-1][-1]=header_all[-1][-1].strip()
                event_id=header_all[-1][-1]
                tmp_pol=[]
                is_header_line=False
            else:
                if line[0:3]=='   ': # Footer line
                    tmp_pick_df=pd.DataFrame(tmp_pol,columns=pick_col_name)

                    # Drops any phase arrivals that do not have a P first motion
                    tmp_pick_df=tmp_pick_df.drop(tmp_pick_df.loc[tmp_pick_df['p_first_motion']==' '].index)
                    # Adds event_id to the dataframe
                    tmp_pick_df['event_id']=str(event_id)

                    pick_all.append(tmp_pick_df)
                    is_header_line=True
                else: # Phase arrival line
                    tmp_pol.append([line[i:j] for i,j in zip(pick_col_ind[:,0], pick_col_ind[:,1])])

    cat_df=pd.DataFrame(header_all,columns=col_name)
    cat_df=cat_df.astype({  'year':'int32','month':'int32',
                            'day':'int32','second':'float',
                            'lat_deg':'int32','lat_min':'int32',
                            'lon_deg':'int32','lon_min':'int32',
                            'depth':'int32','x_error_km':'int32',
                            'z_error_km':'int32','mag_pref':'int32'})

    cat_df.loc[cat_df['year']>50,'year']+=1900
    cat_df.loc[cat_df['year']<50,'year']+=2000
    cat_df.loc[cat_df['hour']=='  ','hour']=0
    cat_df.loc[cat_df['minute']=='  ','minute']=0
    cat_df['hour']=cat_df['hour'].astype('int32')
    cat_df['minute']=cat_df['minute'].astype('int32')
    cat_df['second']/=100
    cat_df['origin_DateTime']=pd.to_datetime(cat_df.loc[:,['year','month','day','hour','minute','second']])
    cat_df['origin_lat']=cat_df['lat_deg']+cat_df['lat_min']/6000
    cat_df.loc[cat_df['latNS']!=' ','origin_lat']*=-1
    cat_df['origin_lon']=cat_df['lon_deg']+cat_df['lon_min']/6000
    cat_df.loc[cat_df['lonEW']==' ','origin_lon']*=-1
    cat_df['origin_depth_km']=cat_df['depth']/100
    cat_df['horz_uncert_km']=cat_df['x_error_km']/100
    cat_df['vert_uncert_km']=cat_df['z_error_km']/100
    cat_df['event_mag']=cat_df['mag_pref']/10

    cat_df=cat_df.loc[:,['origin_DateTime', 'origin_lat', 'origin_lon', 'origin_depth_km',
                        'horz_uncert_km', 'vert_uncert_km', 'event_mag', 'event_id']]

    pick_df=pd.concat(pick_all).reset_index(drop=True)
    pick_df['takeoff_uncertainty']=pick_df['takeoff_uncertainty'].replace('   ',0)
    pick_df['azimuth_uncertainty']=pick_df['azimuth_uncertainty'].replace('   ',0)
    pick_df=pick_df.astype({  'p_weight_code':'int32','sr_dist_km':'float','takeoff':'int32',
                            'azimuth':'int32','takeoff_uncertainty':'int32','azimuth_uncertainty':'int32'})
    pick_df['p_polarity']=0.
    pick_df['sr_dist_km']/=10
    pick_df['p_first_motion']=pick_df['p_first_motion'].str.upper()

    pick_df['takeoff']=180-pick_df['takeoff']

    pick_df.loc[pick_df['p_weight_code']==0,'p_polarity']=1.
    pick_df.loc[pick_df['p_weight_code']==1,'p_polarity']=0.5
    pick_df.loc[pick_df['p_weight_code']==2,'p_polarity']=0.2
    pick_df.loc[pick_df['p_weight_code']==3,'p_polarity']=0.1
    pick_df.loc[(pick_df['p_first_motion']=='D') | (pick_df['p_first_motion']=='-'),'p_polarity']*=-1

    pick_df=pick_df.drop(columns=['p_first_motion','p_weight_code'])

    pick_df['station']=pick_df['station'].str.strip()
    pick_df['channel']=pick_df['channel'].str.strip()

    return cat_df,pick_df


def read_hash3_polarity_file(fpfile,p_weight_I=1.0,p_weight_E=0.5,sta_name_length=4):
    '''
    Reads input phase file following the old SCEC format.
    Gives Impulsive (I) manual picks a weight of p_weight_I (default 1.0)
    Gives Emmergent (E) manual picks a weight of p_weight_E (default 0.5)
    sta_name_length: How many characters the station name has been allocated in the file.
    '''

    col_len=[[0,sta_name_length],
            [sta_name_length+1,sta_name_length+3],
            [sta_name_length+5,sta_name_length+8],
            [sta_name_length+9,sta_name_length+10],
            [sta_name_length+11,sta_name_length+12],
            ]

    file1 = open(fpfile, 'r')
    header_list=[]
    catalog_polarity_list=[]
    event_polarity_list=[]
    is_header=True
    for line in file1:
        if is_header:
            origin_year=int(line[0:4])
            origin_month=int(line[4:6])
            origin_day=int(line[6:8])
            origin_hour=int(line[8:10])
            origin_minute=int(line[10:12])
            origin_second=float(line[12:17])

            origin_lat_deg=int(line[17:19])
            origin_lat_NS=line[19]
            origin_lat_min=float(line[20:25])
            origin_lon_deg=int(line[25:28])
            origin_lon_EW=line[28]
            origin_lon_min=float(line[29:34])

            origin_depth_km=float(line[34:39])
            horz_uncert_km=float(line[88:93])
            vert_uncert_km=float(line[94:99])
            event_mag=float(line[139:143])
            mag_type='M?'
            event_id=str(line[149:165]).strip()

            if origin_lat_NS!='S':
                origin_lat_NS='N'
            if origin_lon_EW!='E':
                origin_lon_EW='W'
            origin_lat=dm2dd(origin_lat_deg,origin_lat_min,origin_lat_NS)
            origin_lon=dm2dd(origin_lon_deg,origin_lon_min,origin_lon_EW)
            header_df=pd.DataFrame( data=[[origin_year,origin_month,origin_day,
                                            origin_hour,origin_minute,origin_second,
                                            origin_lat,origin_lon,origin_depth_km,
                                            horz_uncert_km,vert_uncert_km,
                                            event_mag,mag_type,event_id]],
                                    columns=['year','month','day',
                                            'hour','minute','second',
                                            'origin_lat','origin_lon','origin_depth_km',
                                            'horz_uncert_km','vert_uncert_km',
                                            'event_mag','mag_type','event_id'])
            header_list.append(header_df)
            is_header=False
        elif line[0:4]=='    ':
            tmp_df=pd.DataFrame(data=event_polarity_list,columns=['station','network','location','channel','p_polarity','event_id'])
            try:
                tmp_df=pd.DataFrame(data=event_polarity_list,columns=['station','network','location','channel','p_polarity','event_id'])
            except:
                raise ValueError('*** Error reading line:\n{}'.format(line))
            tmp_df['network']=tmp_df['network'].str.strip()
            tmp_df['station']=tmp_df['station'].str.strip()
            tmp_df['location']=tmp_df['location'].str.strip()
            tmp_df['channel']=tmp_df['channel'].str.strip()
            # tmp_df=tmp_df.sort_values(by=['station','network','channel','location']).reset_index(drop=True)
            catalog_polarity_list.append(tmp_df)
            event_polarity_list=[]
            is_header=True
        else:

            station_name=line[col_len[0][0]:col_len[0][1]]
            network_name=line[col_len[1][0]:col_len[1][1]]
            channel_name=line[col_len[2][0]:col_len[2][1]]
            p_onset=line[col_len[3][0]:col_len[3][1]]
            p_polarity=line[col_len[4][0]:col_len[4][1]]
            location_name='--'
            if (p_onset=='I') | (p_onset=='i'):

                p_weight=p_weight_I
            elif (p_onset=='E') | (p_onset=='e'):
                p_weight=p_weight_E
            else:
                print('Unknown p_onset (',p_onset,'). Discarding polarity.')
                p_weight=0

            if (p_polarity=='U') | (p_polarity=='u') | (p_polarity=='+'):
                p_polarity=p_weight
            elif (p_polarity=='D') | (p_polarity=='d') | (p_polarity=='-'):
                p_polarity=-p_weight
            else:
                print('Unknown p_polarity (',p_polarity,'). Discarding polarity.')
                p_weight=0

            if p_polarity!=0:
                event_polarity_list.append([station_name,network_name,location_name,channel_name,p_polarity,event_id])
    file1.close()

    cat_df=pd.concat(header_list).reset_index(drop=True)
    origin_DateTime=pd.to_datetime(cat_df.loc[:,['year','month','day','hour','minute','second']])
    cat_df.insert(0,'origin_DateTime',origin_DateTime)
    cat_df=cat_df.drop(columns=['year','month','day','hour','minute','second','mag_type'])

    tmp_pol_df=pd.concat(catalog_polarity_list).reset_index(drop=True)
    return cat_df,tmp_pol_df


def read_hash4_polarity_file(fpfile,p_weight_I=1.0,p_weight_E=0.5):
    '''
    Reads a polarity file following the HASH driver 4 format.
    '''
    col_ind_name=np.asarray([
        [0,4,'year'],
        [4,6,'month'],
        [6,8,'day'],
        [8,10,'hour'],
        [10,12,'minute'],
        [12,16,'second'],
        [16,18,'lat_deg'],
        [18,19,'latNS'],
        [19,23,'lat_min'],
        [23,26,'lon_deg'],
        [26,27,'lonEW'],
        [27,31,'lon_min'],
        [31,36,'depth'],
        [130,146,'event_id'],
        [147,150,'mag_pref']
    ])
    col_ind=col_ind_name[:,0:2].astype(int)
    col_name=col_ind_name[:,2]

    pick_col_ind_name=np.asarray([
        [0,4,'station'],
        [5,7,'network'],
        [9,12,'channel'],
        [13,14,'p_onset'],
        [15,16,'p_first_motion'],
        [16,17,'p_weight_code'],
    ])
    pick_col_ind=pick_col_ind_name[:,0:2].astype(int)
    pick_col_name=pick_col_ind_name[:,2]

    header_all=[]
    pick_all=[]
    with open(fpfile, 'r') as file1:
        is_header_line=True
        for line in file1:
            if line[14].upper() in ['P','S','+','-',' ']: # Phase arrival line
                tmp_pol.append([line[i:j] for i,j in zip(pick_col_ind[:,0], pick_col_ind[:,1])])
            else:
                if len(header_all)>0:
                    tmp_pick_df=pd.DataFrame(tmp_pol,columns=pick_col_name)
                    # Drops any phase arrivals that do not have a P first motion
                    tmp_pick_df=tmp_pick_df.drop(tmp_pick_df.loc[tmp_pick_df['p_first_motion']==' '].index)
                    # Adds event_id to the dataframe
                    tmp_pick_df['event_id']=str(event_id)
                    pick_all.append(tmp_pick_df)

                header_all.append([line[i:j] for i,j in zip(col_ind[:,0], col_ind[:,1])])
                header_all[-1][-2]=header_all[-1][-2].strip()
                event_id=header_all[-1][-2]
                tmp_pol=[]
                is_header_line=False
        if len(tmp_pol)>0:
            tmp_pick_df=pd.DataFrame(tmp_pol,columns=pick_col_name)
            # Drops any phase arrivals that do not have a P first motion
            tmp_pick_df=tmp_pick_df.drop(tmp_pick_df.loc[tmp_pick_df['p_first_motion']==' '].index)
            # Adds event_id to the dataframe
            tmp_pick_df['event_id']=event_id
            pick_all.append(tmp_pick_df)

    cat_df=pd.DataFrame(header_all,columns=col_name)
    cat_df=cat_df.astype({  'year':'int32','month':'int32',
                            'day':'int32','second':'float',
                            'lat_deg':'int32','lat_min':'int32',
                            'lon_deg':'int32','lon_min':'int32',
                            'depth':'int32','mag_pref':'int32'})

    cat_df['second']/=100
    cat_df['origin_DateTime']=pd.to_datetime(cat_df.loc[:,['year','month','day','hour','minute','second']])
    cat_df['origin_lat']=cat_df['lat_deg']+cat_df['lat_min']/6000
    cat_df.loc[cat_df['latNS']!=' ','origin_lat']*=-1
    cat_df['origin_lon']=cat_df['lon_deg']+cat_df['lon_min']/6000
    cat_df.loc[cat_df['lonEW']==' ','origin_lon']*=-1
    cat_df['origin_depth_km']=cat_df['depth']/100
    cat_df['event_mag']=cat_df['mag_pref']/100

    cat_df['horz_uncert_km']=0
    cat_df['vert_uncert_km']=0
    cat_df=cat_df.loc[:,['origin_DateTime', 'origin_lat', 'origin_lon', 'origin_depth_km',
                        'horz_uncert_km', 'vert_uncert_km', 'event_mag', 'event_id']]

    pick_df=pd.concat(pick_all).reset_index(drop=True)
    pick_df['p_onset']=pick_df['p_onset'].str.upper()
    pick_df['p_first_motion']=pick_df['p_first_motion'].str.upper()
    pick_df['p_first_motion']=pick_df['p_first_motion'].replace('+','U')
    pick_df['p_first_motion']=pick_df['p_first_motion'].replace('-','D')

    pick_df['p_polarity']=0.
    pick_df.loc[pick_df['p_onset']=='I','p_polarity']=p_weight_I
    pick_df.loc[pick_df['p_onset']=='E','p_polarity']=p_weight_E
    pick_df.loc[(pick_df['p_first_motion']=='D'),'p_polarity']*=-1

    pick_df['station']=pick_df['station'].str.strip()
    pick_df['network']=pick_df['network'].str.strip()
    pick_df['channel']=pick_df['channel'].str.strip()

    pick_df=pick_df.drop(columns=['p_first_motion','p_onset'])

    return cat_df,pick_df
