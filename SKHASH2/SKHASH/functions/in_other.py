'''
Various functions used by SKHASH that do not fall into one of the other file categories.
'''

# External libraries
import numpy as np
import pandas as pd


def read_control_file(control_filepath,p_dict):
    '''
    Reads the user's control file.
    If the control file does not contain a variable, the default value in p_dict is used.
    If the control file contains an unexpected variable, it is ignored after printing a warning message.
    This file should follow the format:
        $var_name
        var_value
        ...
    Where "$var_name" is the name of the variable with a "$" as the first non-whitespace character,
    and "var_value" is the value of the variable on the following line.
    Velocity models ($vmodel_paths) can be specified on multiple lines if desired, e.g.:
        $vmodel_paths
        path_to/vmodel1.txt
        path_to/vmodel2.txt

    User text on all other lines will be ignored, provided the first non-whitespace character on the line is not a "$".
    User text will also be ignored following a space (" ") after "$var_name"
    '''
    # Reads lines from the control file
    with open(control_filepath) as f:
        lines=f.read().splitlines()
        lines=[line.strip() for line in lines]

    # Identifies the variable names and their line numbers
    var_name=[]
    var_line=[]
    value_line=[]
    for line_x,line in enumerate(lines):
        if line: # If not empty
            if line[0]=='$':
                tmp_var_name=line.split()[0][1:]
                tmp_var_name=tmp_var_name.split('#')[0] # In case the user didn't put a space after the variable name.
                if tmp_var_name in p_dict.keys():
                    var_name.append(tmp_var_name)
                    var_line.append(line_x)
                else:
                    raise ValueError('Unknown variable name ({}) in the control file on line {}:\n\t{}'.format(tmp_var_name,line_x,line))
            elif line[0]!='#':
                value_line.append(line_x)

    # Ensures there are not duplicate variables
    var_set=set()
    duplicate_var=[x for x in var_name if x in var_set or var_set.add(x)]
    if duplicate_var:
        raise ValueError('Duplicate variable name in control file: \"{}\"'.format(duplicate_var[0]))

    # Handles velocity model paths that are entered on multiple lines
    if 'vmodel_paths' in var_name:
        tmp_ind=var_name.index('vmodel_paths')
        var_tmp_ind=[x for x in value_line if (x > var_line[tmp_ind])]
        if len(var_tmp_ind)>1:
            if (tmp_ind+1)<len(var_line):
                var_tmp_ind=[x for x in var_tmp_ind if (x < var_line[tmp_ind+1])]
            vmodel_paths_join=' '.join([lines[x] for x in var_tmp_ind])

            value_line=[x for x in value_line if x not in var_tmp_ind[1:]]
            lines[var_tmp_ind[0]]=vmodel_paths_join

    # Ensures that every value line corresponds with a variable
    corr_bool=np.isin( (np.asarray(value_line)-1),var_line)
    if not(np.all(corr_bool)):
        prob_ind=np.where(corr_bool==False)[0][0]
        raise ValueError('Unable to associate a value (\"{}\") on line {} with a variable in the control file. If this is a variable, start the line with a \'$\'. If this is intended to be ignored, start the line with a \'#\'.'.format(lines[value_line[prob_ind]],value_line[prob_ind]))

    # Determines the variable values that follow the variable name lines
    var_value=[]
    for line_x in var_line:
        tmp_var_value=lines[line_x+1].strip()

        if tmp_var_value[0]=='#':
            raise ValueError('Value (\"{}\") on line {} in the control file starts with a \'#\'.'.format(tmp_var_value,line_x+1) )
        if tmp_var_value.isdigit(): # if int
            tmp_var_value=int(tmp_var_value)
        elif tmp_var_value.replace('.','',1).isdigit(): # if float
            tmp_var_value=float(tmp_var_value)
        elif (tmp_var_value=='True'):
            tmp_var_value=True
        elif (tmp_var_value=='False'): # if bool
            tmp_var_value=False
        elif (' ' in tmp_var_value): # if list
            tmp_var_value=tmp_var_value.split()
        var_value.append(tmp_var_value)

    # Assigns the variable values to the dictionary
    for name,value in zip(var_name,var_value):
        p_dict[name]=value

    # Some basic variable type fixing
    p_dict['look_dep']=[int(x) for x in p_dict['look_dep']]
    p_dict['look_del']=[int(x) for x in p_dict['look_del']]
    if type(p_dict['vmodel_paths'])==str:
        p_dict['vmodel_paths']=[p_dict['vmodel_paths']]

    return p_dict


def read_simulps(simulpsfile,pol_df):
    '''
    Appends the takeoffs & azimuths from a SIMULPS 3D (Evans et al., 1994) file
    to the polarity dataframe.
    '''
    duplicated_sta_eventid_df=pol_df[pol_df.duplicated(subset=['station','event_id'])]
    if not(duplicated_sta_eventid_df.empty):
        raise ValueError('When using a SIMULPS file (simulpsfile: {}), the station codes for each event in the polarity files must be unique. Example polarity problems:\n{}'.format(simulpsfile,duplicated_sta_eventid_df.head(10)))

    simulps_df=read_simulps_file(simulpsfile)
    # pol_df=pol_df.merge(simulps_df,on=['station','event_id'])

    duplicated_simulps_df=simulps_df[simulps_df.duplicated(subset=['station','event_id'])]
    if not(duplicated_simulps_df.empty):
        print('*WARNING: The SIMULPS file (simulpsfile: {}) contains duplicate station codes for an individual event. Only the first station record for each event will be considered.'.format(simulpsfile))
        simulps_df=simulps_df.drop(duplicated_simulps_df.index).reset_index(drop=True)

    # Updates the source-receiver distances, takeoff, and azimuths. Information in SIMULPS file will overwrite any previous distances/angles.
    pol_df.update(pol_df.drop(columns=['sr_dist_km', 'azimuth', 'takeoff', 'azimuth_uncertainty', 'takeoff_uncertainty']).merge(simulps_df, how='left', on=['station','event_id']))

    return pol_df


def read_simulps_file(simulpsfile):
    '''
    Reads a SIMULPS 3D (Evans et al., 1994) file. See HASH v1.2 manual for the required formatting.
    '''
    col_ind_name=np.asarray([
        [1,3,'year'],
        [3,5,'month'],
        [5,7,'day'],
        [8,10,'hour'],
        [10,12,'minute'],
        [13,18,'second'],
        [19,21,'lat_deg'],
        [21,22,'latNS'],
        [22,27,'lat_min'],
        [28,31,'lon_deg'],
        [31,31,'lonEW'],
        [32,37,'lon_min'],
        [38,44,'depth'],
        [47,50,'mag_pref'],
        [55,63,'event_id'],
    ])
    col_ind=col_ind_name[:,0:2].astype(int)

    pick_col_ind_name=np.asarray([
        [1,5,'station'],
        [6,11,'sr_dist_km'],
        [11,15,'azimuth'],
        [15,19,'takeoff'],
        [65,69,'azimuth_uncertainty'],
        [69,72,'takeoff_uncertainty'],
    ])
    pick_col_ind=pick_col_ind_name[:,0:2].astype(int)
    pick_col_name=pick_col_ind_name[:,2]

    with open(simulpsfile, 'r') as file1:
        lines = file1.read().splitlines()
    header_line_nums=[]
    station_startline_nums=[]
    station_endline_nums=[]
    station_line_flag=False
    for line_x in range(len(lines)):
        if lines[line_x][:16]=='  DATE    ORIGIN':
            header_line_nums.append(line_x)
        elif lines[line_x][:11]=='  STN  DIST':
            station_startline_nums.append(line_x+1)
            station_line_flag=True
        elif station_line_flag & (lines[line_x].strip()==''):
            station_endline_nums.append(line_x-1)
            station_line_flag=False
    if station_line_flag:
        station_endline_nums.append(line_x-1)

    if len(station_startline_nums)!=len(station_endline_nums):
        raise ValueError('Error reading simulpsfile ({}).'.format(simulpsfile))

    header_all=[]
    event_id=[]
    for header_x in range(len(header_line_nums)):
        line=lines[header_line_nums[header_x]+1]
        header_all.append([line[i:j] for i,j in zip(col_ind[:,0], col_ind[:,1])])
        event_id.append(header_all[-1][-1].strip())

    simulps_all=[]
    for event_x in range(len(station_startline_nums)):
        tmp_event=[]
        for line_x in range(station_startline_nums[event_x],station_endline_nums[event_x]+1):
            tmp_event.append([lines[line_x][i:j] for i,j in zip(pick_col_ind[:,0], pick_col_ind[:,1])])
        tmp_pick_df=pd.DataFrame(tmp_event,columns=pick_col_name)
        tmp_pick_df['event_id']=event_id[event_x]
        simulps_all.append(tmp_pick_df)
    simulps_df=pd.concat(simulps_all).reset_index(drop=True)

    simulps_df=simulps_df.astype({'sr_dist_km':float,'azimuth':int,'takeoff':int,'azimuth_uncertainty':float,'takeoff_uncertainty':float,'event_id':str})
    simulps_df['station']=simulps_df['station'].str.strip()
    simulps_df['takeoff']=180-simulps_df['takeoff']
    return simulps_df


def read_catalog_file(p_dict):
    '''
    Reads catalog file
    '''
    cat_df=pd.read_csv(p_dict['catfile'],skipinitialspace=True,comment='#')
    if not({'time','event_id','latitude','longitude','depth'}.issubset(cat_df.columns)):
        raise ValueError(('The catalog file (catfile: {}) must contain the following column names: [time, event_id, latitude, longitude, depth]\n'+
                        'Only the following columns were provided:\n{}').format(p_dict['catfile'],cat_df.columns.values))
    cat_df=cat_df.rename(columns={'time':'origin_DateTime','latitude':'origin_lat','longitude':'origin_lon','depth':'origin_depth_km'})

    if not('vert_uncert_km' in cat_df.columns):
        cat_df['vert_uncert_km']=p_dict['default_vert_uncert_km']
    if not('horz_uncert_km' in cat_df.columns):
        cat_df['horz_uncert_km']=p_dict['default_horz_uncert_km']

    if 'origin_DateTime' in cat_df.columns:
        cat_df['origin_DateTime']=pd.to_datetime(cat_df['origin_DateTime'])
        cat_df=cat_df[['origin_DateTime','origin_lat','origin_lon','origin_depth_km','horz_uncert_km','vert_uncert_km','event_id']]
    else:
        cat_df=cat_df[['origin_lat','origin_lon','origin_depth_km','horz_uncert_km','vert_uncert_km','event_id']]

    # Ensures column data formats are the expected type
    cat_df['origin_lat']=cat_df['origin_lat'].astype(float)
    cat_df['origin_lon']=cat_df['origin_lon'].astype(float)
    cat_df['origin_depth_km']=cat_df['origin_depth_km'].astype(float)
    cat_df['horz_uncert_km']=cat_df['horz_uncert_km'].astype(float)
    cat_df['vert_uncert_km']=cat_df['vert_uncert_km'].astype(float)
    cat_df['event_id']=cat_df['event_id'].astype(str)

    cat_df=cat_df.sort_values(by='event_id').reset_index(drop=True)
    
    return cat_df
