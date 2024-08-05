'''
Python package for earthquake focal mechanism inversions.

Author: Robert J. Skoumal (rskoumal@usgs.gov)

This project contains Python code to compute focal mechanism solutions using
first-motion polarities (traditional, consensus, and/or imputed) and S/P ratios
(traditional and/or consensus).

Significant portions of this code are based on the HASH algorithm created by
Jeanne L. Hardebeck and Peter M. Shearer. If you use this code, please cite the
appropriate papers.
	Skoumal, R.J., Hardebeck, J.L., & Shearer, P.M. (in press). SKHASH: A Python
		package for computing earthquake focal mechanisms. Seismological Research Letters.
	Hardebeck, J.L., & Shearer, P.M. (2002). A new method for determining
		first-motion focal mechanisms. Bulletin of the Seismological Society of America,
		92(6), 2264-2276. https://doi.org/10.1785/0120010200
	Hardebeck, J.L., & Shearer, P.M. (2003). Using S/P amplitude ratios to
		constrain the focal mechanisms of small earthquakes. Bulletin of the Seismological
		Society of America, 93(6), 2434-2444. https://doi.org/10.1785/0120020236
'''

# Standard libraries
import sys
import time
import importlib.util
import multiprocessing
import argparse

# External libraries
import numpy as np
import pandas as pd

# Local libraries
import functions.in_pol as in_pol # Reading polarity inputs
import functions.in_sta as in_sta # Reading station files
import functions.in_sp as in_sp # Reading S/P ratios
import functions.in_qc as in_qc # Quality control of user inputs
import functions.in_other as in_other # Functions for reading other inputs
import functions.fun as fun # Computing mechanisms
import functions.out as out # Output functions
import functions.gridsearch_so as gridsearch_so # For creating fortran gridsearch module
import functions.compute_mech as compute_mech # For computing mechanism

# Superficial version information
version_string='v0.1'
version_date='2024-02-16'

'''
Sets default parameter values. If a value is not provided in the control file (or command line) for a variable, the default value is used.
'''
p_dict={
	'input_format':'skhash', # default input file format
	'input_format_stfile':'', # station list file format
	'input_format_fpfile':'', # traditional P-polarity file format
	'input_format_impfile':'', # imputed P-polarity input
	'input_format_conpfile':'', # consensus P-polarities file format
	'input_format_dlpfile':'', # consensus P-polarities file format
	'input_format_ampfile':'', # amplitude input file format
	'input_format_relampfile':'', # relative S/P ratio file format
	'input_format_plfile':'', # polarity reversal format

	'controlfile':'',
	'catfile':'', # earthquake catalog information filename
	'stfile':'', # station list filename
	'plfile':'', # station polarity reversal filename
	'corfile':'', # station correction filename
	'fpfile':'', # traditional P-polarity input filename
	'impfile':'', # imputed P-polarity input filename
	'conpfile':'', # consensus P-polarities from relative measurements input filename
	'dlpfile':'', # consensus P-polarities from relative measurements input filename
	'ampfile':'', # amplitude input filename
	'relampfile':'', # relative S/P ratios input filename
	'simulpsfile':'', # filename for SIMULPS azimuth and takeoff angles
	'outfile1':'default_out1.txt', # focal mechanisms output filename
	'outfile2':'', # acceptable plane output filename
	'outfile_pol_agree':'', # record of polarity (dis)agreeement output filename
	'outfile_sp_agree':'', # record of S/P difference output filename
	'outfile_pol_info':'', # record of all polarities considered in the mechanisms
	'outfolder_plots':'', # Folder where simple focal mechanism plots will be created (outfolder_plots/event_id.png). To ignore, leave blank.

	'npolmin':8, # mininum number of polarity data (e.g., 8)
	'nmc':30, # number of trials (e.g., 30)
	'maxout':500, # maxout for focal mech. output (e.g., 500)
	'ratmin':3.0, # minimum allowed signal-to-noise ratio
	'badfrac':0.1, # fraction polarities assumed bad
	'badmin':2.0, # minimum number of polarities per event assumed bad. If badmin > [(sum of the polarity weights) * badfrac], badmin replaces badfrac
	'qbadfrac':0.3, # assumed noise in amplitude ratios, log10 (e.g. 0.3 for a factor of 2)
	'qbadmin':2.0, # minimum S/P amplitude misfit allowed per event. If qbadmin > [(sum of the number of S/P measurements) * qbadfrac], qbadmin replaces qbadfrac
	'min_polarity_weight':0.1, # Any polarities with an abs(weight) < min_polarity_weight will be ignored
	'delmin':0.0, # minimum allowed source-station distance in km. If delmin<=0, uses all results >= than look_del[0] distance.
	'delmax':120.0, # maximum allowed source-station distance in km. If delmax<=0, uses all results <= than look_del[1] distance.
	'azmax':10.0, # maximum allowed source-station azimuth uncertainty in degrees. To use all results regardless of azimuthal uncertainty, set azmax=0.
	'pmax':0.0, # maximum allowed takeoff ("plungal") uncertainty in degrees. To use all results regardless of takeoff uncertainty, set pmax=0.
	'cangle':45.0, # angle for computing mechanisms probability
	'prob_max':0.2, # probability threshold for multiples (e.g., 0.2)
	'max_agap':90.0, # maximum azimuthal gap
	'max_pgap':60.0, # maximum "plungal" gap

	'vmodel_paths':[], # list of paths to velocity model
	'write_lookup_table':False, # The lookup table will be saved so that it can be used in future runs. A .npy file extension will be added to the vmodel_paths.
	'recompute_lookup_table':True, # If a velocity model with an identical filename already exists, it will remake it.

	'require_temporal_match':False, # Requires the reported earthquake to occur between the station metadata start-end time.
	'require_network_match':True, # Requires the reported pick network to match the metadata
	'require_station_match':True, # Requires the reported pick station to match the metadata
	'require_channel_match':True, # Requires the reported pick channel to match the metadata
	'require_location_match':True, # Requires the reported pick location to match the metadata
	'allow_duplicate_stations':False, # Allows an event_id to have multiple polarities and S/P ratios from the same receiver
	'iterative_avg':False, # If True, it will iteratively compute the average mech by removing the solution furthest from the avg following HASH.

	'min_quality_report':'', # Only mech qualities of this or better will be used to reported. To consider all accepted solutions, leave blank.

	'overwrite_output_file':True, # If True, SKHASH will overwrite any existing output files.
	'perturb_epicentral_location':False, # If True, it will randomly perturb the horizontal location of the earthquakes

	'ignore_missing_metadata':False, # If True, any measurements with missing metadata will be ignored. If False, an error will be raised.

	'use_fortran':False, # If True, it will attempt to load the Fortran subroutine. If it cannot load or is False, the Python routine will be used.
	'num_cpus':1, # Number of cores to run things in parallel. Set to 0 to use all available, or 1 to run in serial.

	'dang':5, # minimum grid spacing (degrees).
	'nx0':101, # maximum source-station distance bins for look-up tables
	'nd0':14, # maximum source depth bins for look-up tables
	'look_dep':[0,39,3], # minimum source depth, maximum, and interval for the lookup table
	'look_del':[0,200,2], # minimum source-station distance, maximum, and interval for the lookup table
	'allow_hypocenters_outside_table':False, # If a hypocenter is outside of the lookup table range, the min/max lookup value is used. If False, an error will be produced.
	'nump':9000, # number of rays traced

	'min_sp':0.0005,# If any user-provided S/P ratios are less than min_sp after applying station corrections, those values are replaced with min_sp. To ignore, use a value <= 0.
	'max_sp':2000.0, # If any user-provided S/P ratios are greater than max_sp after applying station corrections, those values are replaced with max_sp. To ignore, use a value <= 0.
	'min_amp':0.0005, # minimum amplitude [0-1] considered when calculating expected P and S amplitudes. Values < than this are treated as 0. To ignore, use a value <= 0.

	'epicenter_degree_precision':6, # Number of decimal places to round lat, lon epicenter coordinates.
	'output_angle_precision':1, # Number of decimal places to output for angles (e.g., strike, dip, rake, takeoff)
	'output_quality_precision':1, # Number of decimal places to output for mech qualities (e.g., prob, mfrac, mavg, stdr)
	'output_km_distance_precision':3, # Number of decimal places to output for source-receiver distances (reported in km)
	'output_vector_precision':4, # Number of decimal places to output for normal vectors
	'default_vert_uncert_km':1.0, # If the vertical uncertainty for an earthquake is missing or negative, this value is used instead when perturbing the source depth.
	'default_horz_uncert_km':1.0, # If the horizontal uncertainty for an earthquake is missing or negative, this value is used instead when perturbing the source hypocenter.

	'plot_station_names':False, # If True, station names will be included in the beachball plot
	'plot_acceptable_solutions':False, # If True, all acceptable solutions will be plotted in the beachball plot

	# Warning: only change the following values if you know what you're doing! :)
	'compute_takeoff_azimuth':True, # Used to determine if takeoff and source-receiver azimuths need to be computed
	'npick0':15000, # maximum number of picks per event. Used by Fortran gridsearch. If this is changed, the value in gridsearch.f must be changed and recompiled.
	'merge_on':[] # Controls how to merge metadata with observations. List that includes ('network', 'station', 'location', 'channel')
}

# Quality codes should be listed from best quality (most restrictive) to worst quality (least restrictive)
qual_criteria_dict={
				'qual_letter':np.asarray([ 'A',  'B',  'C',     'D']),
				'probs':np.asarray(      [ 0.8,  0.6,  0.5,      0 ]),
				'var_avg':np.asarray(    [  25,   35,   45, np.inf ]),
				'mfrac':np.asarray(      [0.15, 0.20, 0.30, np.inf ]),
				'stdr':np.asarray(       [ 0.5,  0.4,  0.3,      0 ]),
				}

'''
Reads command line arguments.
Note that if you also use a control file, the control file variable will overwrite command line arguments.
'''
parser = argparse.ArgumentParser()
parser.add_argument('controlfile',nargs='?')
for p_dict_var in p_dict.keys():
	parser.add_argument('--'+p_dict_var)
args = vars(parser.parse_args())
for p_dict_var in p_dict.keys():
	if args[p_dict_var] is not None:
		dtype=type(p_dict[p_dict_var])
		p_dict[p_dict_var]=dtype(args[p_dict_var])

if __name__ == "__main__":
	print('========================\nSKHASH {} ({})\n========================'.format(version_string,version_date))
	total_runtime_start=time.time()

	'''
	Checks Python version.
	'''
	if not((sys.version_info.major==3) & (sys.version_info.minor>=8)):
		print('***WARNING: SKHASH has only been tested on Python 3.8+. Refer to the manual for supported versions.')

	'''
	Reads the user-created control file, overwritting the default and commonand line parameter values.
	'''
	if p_dict['controlfile']:
		print('Control file: {}'.format(p_dict['controlfile']))
		p_dict=in_other.read_control_file(p_dict['controlfile'],p_dict)

	'''
	Loads/creates the fortran extension module
	'''
	if p_dict['use_fortran']:
		gridsearch_spec = importlib.util.find_spec("functions.gridsearch")
		if (gridsearch_spec is not None):
			import functions.gridsearch as gridsearch
			print('Fortran subroutine module successfully loaded.')
		else:
			status=gridsearch_so.create()
			if status:
				import functions.gridsearch as gridsearch
			else:
				p_dict['use_fortran']=False
				print('Unable to load Fortran subroutine. Falling back to using the Python gridsearch routine.\n')

	'''
	Performs various checks to ensure that the user-defined parameters are appropriate
	'''
	p_dict=in_qc.check_input_params(p_dict,qual_criteria_dict)

	'''
	Reads P-wave first motion polarity file(s)
	'''
	pol_df=[]
	cat_df=[]
	if p_dict['fpfile']:
		tmp_cat_df,tmp_pol_df=in_pol.read_polarity_file(p_dict['fpfile'],p_dict['input_format_fpfile'],p_dict['merge_on'])
		tmp_pol_df['source']='trad_p'
		if len(tmp_pol_df):
			pol_df.append(tmp_pol_df)
		if len(tmp_cat_df):
			cat_df.append(tmp_cat_df)
	if p_dict['impfile']:
		tmp_cat_df,tmp_pol_df=in_pol.read_polarity_file(p_dict['impfile'],p_dict['input_format_impfile'],p_dict['merge_on'])
		tmp_pol_df['source']='imp_p'
		if len(tmp_pol_df):
			pol_df.append(tmp_pol_df)
		if len(tmp_cat_df):
			cat_df.append(tmp_cat_df)
	if p_dict['conpfile']:
		tmp_cat_df,tmp_pol_df=in_pol.read_polarity_file(p_dict['conpfile'],p_dict['input_format_conpfile'],p_dict['merge_on'])
		tmp_pol_df['source']='con_p'
		if len(tmp_pol_df):
			pol_df.append(tmp_pol_df)
		if len(tmp_cat_df):
			cat_df.append(tmp_cat_df)
	if p_dict['dlpfile']:
		tmp_cat_df,tmp_pol_df=in_pol.read_polarity_file(p_dict['dlpfile'],p_dict['input_format_dlpfile'],p_dict['merge_on'])
		tmp_pol_df['source']='dl_p'
		if len(tmp_pol_df):
			pol_df.append(tmp_pol_df)
		if len(tmp_cat_df):
			cat_df.append(tmp_cat_df)			
	if len(pol_df):
		pol_df=pd.concat(pol_df).reset_index(drop=True)
	else:
		pol_df=pd.DataFrame()
		print('*WARNING: No polarity measurements provided.')
	if len(cat_df):
		cat_df=pd.concat(cat_df).drop_duplicates(subset='event_id').reset_index(drop=True)

	# Looks for duplicate measurements
	if not(p_dict['allow_duplicate_stations']):
		if len(pol_df):
			duplicate_polarity_index=pol_df[pol_df.duplicated()].index
			if not(duplicate_polarity_index.empty):
				raise ValueError('Duplicate information for {} measurements. Example issues:\n{}'.
					format(len(duplicate_polarity_index),pol_df.loc[duplicate_polarity_index,:].head(10).to_string()))

	# Adds SIMULPS source-receiver distances, takeoff, and azimuths
	if p_dict['simulpsfile']:
		pol_df=in_other.read_simulps(p_dict['simulpsfile'],pol_df)

	'''
	Reads S/P amplitude file, calculates S/P ratios, discards S/P ratios below the noise threshold, and applies station corrections.
	'''
	if p_dict['ampfile']:
		spamp_df,p_dict['ampfile']=in_sp.read_amp_corr_files(p_dict['ampfile'],p_dict)
		spamp_df['source']='trad_sp'

		# Concats polarity and S/P data into a single dataframe
		pol_df=pd.concat([pol_df,spamp_df]).reset_index(drop=True)

	if p_dict['relampfile']:
		spamp_df,p_dict=in_sp.read_amp_corr_files(p_dict['relampfile'],p_dict)
		spamp_df['source']='rel_sp'

		# Concats polarity and consensus S/P data into a single dataframe
		pol_df=pd.concat([pol_df,spamp_df]).reset_index(drop=True)

	'''
	Reads earthquake catalog
	'''
	if p_dict['catfile']:
		cat_df=in_other.read_catalog_file(p_dict)

	# Adds event times, locations, and uncertainties to polarity information
	if len(cat_df):
		if not('event_id2' in pol_df.columns):
			pol_df=pol_df.merge(cat_df,on='event_id',how='left')
		else: # composite mechs
			pol_df=pol_df.merge(cat_df.rename(columns={'event_id':'event_id2'}),on='event_id2',how='left')

	'''
	Quality controls the datasets
	'''
	# Quality control of catalog dataset
	cat_df=in_qc.check_cat_df(cat_df,pol_df,p_dict)

	# Quality control of polarity dataset
	pol_df=in_qc.check_pol_df(pol_df,cat_df,p_dict)

	'''
	Checks that earthquake locations will be contained in the model
	'''
	if p_dict['compute_takeoff_azimuth']:
		pol_df,p_dict=in_qc.check_model_eq(cat_df,pol_df,p_dict)

	'''
	Reads station polarity reversals and applies them to the polarity measurements
	'''
	if p_dict['plfile']:
		pol_reverse_df=in_sta.read_reverse_file(p_dict)
		pol_df=in_sta.reverse_polarities(pol_df,pol_reverse_df,p_dict)
	pol_df=pol_df.drop(pol_df.filter(['station']),axis=1)

	'''
	Reads station metadata file and appends the locations to the polarities
	'''
	if p_dict['stfile']:
		station_df=in_sta.read_station_file(pol_df,p_dict)
		pol_df=in_sta.apply_station_locations(pol_df,station_df,p_dict)

	# Drops columns that are no longer needed
	if 'origin_DateTime' in pol_df:
		pol_df=pol_df.drop(columns=['origin_DateTime'])

	'''
	Dropping any measurements with source-receiver distances > delmax
	'''
	pol_df['sr_dist_km']=pol_df['sr_dist_km'].round(p_dict['output_km_distance_precision'])
	if (p_dict['delmax']>0) | (p_dict['delmin']>0):
		drop_flag=((pol_df['sr_dist_km']>p_dict['delmax']) | (pol_df['sr_dist_km']<p_dict['delmin']))
		pol_df=pol_df.drop(pol_df[drop_flag].index).reset_index(drop=True)
		if np.any(drop_flag):
			if p_dict['delmin']>0:
				print('Discarded {} of the {} polarity measurements with source-receiver distances <{} or >{} km'.
					format(np.sum(drop_flag),len(pol_df),p_dict['delmin'],p_dict['delmax']))
			else:
				print('Discarded {} of the {} polarity measurements with source-receiver distances >{} km'.
					format(np.sum(drop_flag),len(pol_df),p_dict['delmax']))

	'''
	Calculates maximum possible source-receiver azimuth variation, dropping any polarities > azmax
	'''
	if p_dict['azmax']>0:
		if not('azimuth_uncertainty' in pol_df.columns):
			sr_unc_flag=pol_df['horz_uncert_km']<=pol_df['sr_dist_km']

			pol_df['azimuth_uncertainty']=(90-np.rad2deg(np.arccos(pol_df['horz_uncert_km']/pol_df['sr_dist_km'].values)))*2

			# If the horizontal uncertainty is greater than the source-receiver distance, any azimuth is possible.
			tmp_index=pol_df[~sr_unc_flag].index
			pol_df.loc[tmp_index,'azimuth_uncertainty']=360

		# Drop polarities that have a source-receiver azimuth variation > azmax
		tmp_index=pol_df[pol_df['azimuth_uncertainty'] > p_dict['azmax']].index
		if not(tmp_index.empty):
			pol_df.drop(tmp_index, inplace=True)
			print('Discarded {} polarity measurement(s) with source-receiver azimuth uncertainties >{} deg'.format(len(tmp_index),p_dict['azmax']))

		pol_df['azimuth_uncertainty']=pol_df['azimuth_uncertainty'].round(p_dict['output_angle_precision'])
		pol_df=pol_df.reset_index(drop=True)

	'''
	If using precomputed takeoff/azimuths, discards polarities with takeoffs with large uncertainties.
	If computing takeoff/azimuths, this is handled later in compute_mech().
	'''
	# Discards any measurements with a takeoff uncertainty > pmax
	if not(p_dict['stfile']):
		if p_dict['pmax']>0:
			if 'takeoff_uncertainty' in pol_df.columns:
				drop_index=pol_df[(pol_df['takeoff_uncertainty']>p_dict['pmax'])].index
				if len(drop_index)>0:
					print('Discarded {} polarity measurement(s) with takeoff uncertainties >{} deg'.format(len(drop_index),p_dict['pmax']))
					pol_df=pol_df.drop(drop_index).reset_index(drop=True)

	'''
	Discards earthquakes with fewer than npolmin P-polarities
	'''
	if p_dict['npolmin']>0:
		tmp_df=pol_df.groupby(by='event_id').agg({'p_polarity':'count'})
		discard_event_id=tmp_df[tmp_df['p_polarity']<p_dict['npolmin']].index
		if not(discard_event_id.empty):
			pol_df=pol_df.drop(pol_df[pol_df['event_id'].isin(discard_event_id)].index).reset_index(drop=True)
			print('Discarded {} earthquakes with fewer than {} polarities.'.format(len(discard_event_id),p_dict['npolmin']))

	'''
	Removes cataloged earthquakes that have no selected measurements
	'''
	if len(cat_df):
		cat_consider_flag=~(cat_df['event_id'].isin(pol_df['event_id']))
		if cat_consider_flag.any():
			cat_df=cat_df.drop(cat_df[cat_consider_flag].index).reset_index(drop=True)

	'''
	Reads the velocity model files and creates (or loads) the lookup tables.
	'''
	if p_dict['compute_takeoff_azimuth']:
		lookup_dict=fun.create_lookup_table(p_dict)
	else:
		lookup_dict={'deptab':[],'delttab':[],'table':[]}

	'''
	Sets up array with direction cosines for all coordinate transformations
	'''
	if p_dict['use_fortran']:
		dir_cos_dict={}
	else:
		dir_cos_dict=fun.dir_cos_setup(p_dict)

	'''
	Groups polarities and S/P ratios by event_id
	'''
	group_pol_df=pol_df.groupby(by='event_id')
	event_ids=list(group_pol_df.groups.keys())

	if p_dict['outfile_pol_agree']:
		pol_df['pol_agreement']=0
	if p_dict['outfile_sp_agree']:
		pol_df['sp_diff']=-999.

	# Creates outfile1/outfile2 and adds header lines to them.
	if p_dict['outfile1']:
		out.create_outfile1(p_dict['outfile1'],cat_df,pol_df)
	if p_dict['outfile2']:
		out.create_outfile2(p_dict['outfile2'])

	if not(event_ids):
		print('No mechanisms to compute. Exiting.')
		quit()

	'''
	Computes the focal mechanisms
	'''
	mech_runtime_start=time.time()
	num_events=len(event_ids)
	if p_dict['num_cpus']==1: # Run in serial
		print('Computing mechanisms in serial...')
		for event_x,event_id in enumerate(event_ids):
			mech_dict=compute_mech.compute_mech(event_x,num_events,event_id,group_pol_df.get_group(event_id),
												p_dict,lookup_dict,qual_criteria_dict,cat_df,dir_cos_dict)
			if mech_dict['mech_qual']:
				pol_df.loc[mech_dict['event_index'],'mech_quality']=mech_dict['mech_qual']
				if p_dict['outfile_pol_agree']:
					pol_df.loc[mech_dict['event_index'],'pol_agreement']=mech_dict['pol_agreement_out']
				if p_dict['outfile_sp_agree']:
					pol_df.loc[mech_dict['event_index'],'sp_diff']=mech_dict['sp_diff_out']
				if (p_dict['stfile']) and (p_dict['outfile_pol_info']):
					pol_df.loc[mech_dict['event_index'],'takeoff']=mech_dict['takeoff']
					pol_df.loc[mech_dict['event_index'],'azimuth']=mech_dict['sr_az']
					pol_df.loc[mech_dict['event_index'],'takeoff_uncertainty']=mech_dict['takeoff_uncertainty']
					pol_df.loc[mech_dict['event_index'],'azimuth_uncertainty']=mech_dict['azimuth_uncertainty']
	else: # Run in parallel
		print('Computing mechanisms in parallel...')
		pool=multiprocessing.Pool(processes=p_dict['num_cpus'])

		async_results=[]
		for event_x,event_id in enumerate(event_ids):
			try:
				event_pol_df=group_pol_df.get_group(event_id)
			except:
				print('Error getting parallel result for event_id: {}'.format(event_id))
				continue
			async_results.append(pool.apply_async(compute_mech.compute_mech,
						args=(event_x,num_events,event_id,event_pol_df,p_dict,lookup_dict,qual_criteria_dict,cat_df,dir_cos_dict)))
		pool.close()
		pool.join()

		if any([p_dict['outfile_pol_agree'],p_dict['outfile_sp_agree'],p_dict['outfile_pol_info']]):
			for result in async_results:
				mech_dict=result.get()
				if mech_dict['mech_qual']:
					pol_df.loc[mech_dict['event_index'],'mech_quality']=mech_dict['mech_qual']
					if p_dict['outfile_pol_agree']:
						pol_df.loc[mech_dict['event_index'],'pol_agreement']=mech_dict['pol_agreement_out']
					if p_dict['outfile_sp_agree']:
						pol_df.loc[mech_dict['event_index'],'sp_diff']=mech_dict['sp_diff_out']
					if (p_dict['stfile']) and (p_dict['outfile_pol_info']):
						pol_df.loc[mech_dict['event_index'],'takeoff']=mech_dict['takeoff']
						pol_df.loc[mech_dict['event_index'],'azimuth']=mech_dict['sr_az']
						pol_df.loc[mech_dict['event_index'],'takeoff_uncertainty']=mech_dict['takeoff_uncertainty']
						pol_df.loc[mech_dict['event_index'],'azimuth_uncertainty']=mech_dict['azimuth_uncertainty']

	print('Mech computation runtime: {:.2f} sec'.format(time.time()-mech_runtime_start), flush=True)

	'''
	Determines the polarity agreements at the different stations and writes it to file
	'''
	if p_dict['outfile_pol_agree']:
		out.pol_agree(pol_df,p_dict)

	'''
	Determines the S/P agreements at the different stations and writes it to file
	'''
	if p_dict['outfile_sp_agree']:
		out.sp_agree(pol_df,p_dict)

	'''
	Creates output file with lots of information about the polarity measurements used to compute the focal mechanisms
	'''
	if p_dict['outfile_pol_info']:
		out.pol_info(pol_df,p_dict)

	print('Total runtime: {:.2f} sec'.format(time.time()-total_runtime_start), flush=True)
