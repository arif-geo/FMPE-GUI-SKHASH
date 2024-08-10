
# Standard libraries
import time

# External libraries
import numpy as np

import functions.gridsearch_so as gridsearch_so # For creating fortran gridsearch module
import functions.out as out # Output functions
import functions.fun as fun # Computing mechanisms


def compute_mech(event_x,num_events,event_id,event_pol_df,p_dict,lookup_dict,qual_criteria_dict,cat_df,dir_cos_dict):
    '''
    Computes focal mechanisms.
    Input:
        event_x: The event number (cosmetic)
        num_events: Total number of events to compute mechanisms (cosmetic)
        event_id: event id string for the event
        event_pol_df: polarity dataframe
        p_dict: Parameter values created in SKHASH.py, dictionary
        lookup_dict: dictionary with lookup variables, produced by create_lookup_table()
        qual_criteria_dict: dictionary of quality criteria, created in SKHASH.py
        cat_df: catalog dataframe
        dir_cos_dict: dictionary of coordinate transformation variables, created by dir_cos_setup()
    Output:
        mech_dict: dictionary of mechanism solutions
    '''
    event_runtime_start = time.time()

    mech_dict={'event_index':-1,'pol_agreement_out':[],
            'sp_diff_out':[],'takeoff':-1.,'sr_az':-1.,
            'takeoff_uncertainty':-1.,'azimuth_uncertainty':-1.,
            'mech_qual':''}

    if p_dict['stfile']: # Perturb earthquake locations and determine azimuth and takeoff angles
        perturbed_origin_depth_km,sr_dist_km,sr_azimuth=fun.perturb_eq_locations(event_pol_df,p_dict['look_dep'],p_dict['perturb_epicentral_location'],nmc=p_dict['nmc'])
        takeoff=fun.lookup_takeoff(lookup_dict['table'],perturbed_origin_depth_km,sr_dist_km,p_dict['look_dep'],p_dict['look_del'],lookup_dict['deptab'],lookup_dict['delttab'],num_velocity_models=len(p_dict['vmodel_paths']))

        # Discards any measurements with a takeoff uncertainty > pmax
        if p_dict['pmax']>0:
            calc_takeoff_uncertainty=np.std(takeoff,axis=1)
            rm_ind=np.where(calc_takeoff_uncertainty>p_dict['pmax'])[0]

            event_pol_df=event_pol_df.drop(event_pol_df.index[rm_ind])
            sr_azimuth=np.delete(sr_azimuth,rm_ind,axis=0)
            takeoff=np.delete(takeoff,rm_ind,axis=0)
            print('Dropped {} measurements for high takeoff uncertainties'.format(len(rm_ind)))
    else: # Perturb predetermined azimuth and takeoff angles
        sr_azimuth,takeoff=fun.perturb_azimuth_takeoff(event_pol_df,p_dict['nmc'])

    # Calculates the maximum azimuthal and takeoff angle gaps, skipping the event if necessary
    max_azimuthal_gap,max_takeoff_gap=fun.determine_max_gap(sr_azimuth[:,0],takeoff[:,0])

    if max_azimuthal_gap>p_dict['max_agap']:
        print('{} / {}\t({})\n'.format(event_x,num_events-1,event_id)+
              '\tMaximum azimuthal gap ({}) > max_agap ({}). Skipping.'.format(max_azimuthal_gap.round(3),p_dict['max_agap']))
        return mech_dict
    if max_takeoff_gap>p_dict['max_pgap']:
        print('{} / {}\t({})\n'.format(event_x,num_events-1,event_id)+
              '\tMaximum takeoff angle gap ({}) > max_pgap ({}). Skipping.'.format(max_takeoff_gap.round(3),p_dict['max_pgap']))
        return mech_dict

    # P-polarity parameters for determining best-fit solutions
    p_pol=event_pol_df['p_polarity'].values
    p_pol[~np.isfinite(p_pol)]=0
    sumpolweight=np.sum(np.abs(p_pol))
    nextra=max([round(sumpolweight*p_dict['badfrac']*0.5),p_dict['badmin']])
    ntotal=max([round(sumpolweight*p_dict['badfrac']),p_dict['badmin']])

    # S/P ratio parameters for determining best-fit solutions
    if p_dict['ampfile'] or p_dict['relampfile']:
        sp_amp=event_pol_df['sp_ratio'].values
        nspr=sum(np.isfinite(sp_amp))

        # additional amplitude misfit allowed above minimum
        qextra=max([nspr*p_dict['qbadfrac']*0.5,p_dict['qbadmin']])

        # total allowed amplitude misfit
        qtotal=max([nspr*p_dict['qbadfrac'],p_dict['qbadmin']])
    else:
        sp_amp=np.empty(len(event_pol_df));sp_amp[:]=np.nan
        qextra=0
        qtotal=0

    # Runs the gridsearch to find potential mech solutions.
    if p_dict['use_fortran']: # Uses the Python C/API to call Fortran subroutine
        import functions.gridsearch as gridsearch
        p_azi_mc,p_the_mc,f_sp_amp,f_p_pol,p_qual=gridsearch_so.prep_subroutine(sr_azimuth,takeoff,p_pol,sp_amp,p_dict['npick0'],p_dict['nmc'])
        nf,strike_all,dip_all,rake_all,faultnorms_all,faultslips_all=gridsearch.focalamp_mc_wt(p_azi_mc,p_the_mc,f_sp_amp,f_p_pol,p_qual,p_dict['nmc'],p_dict['dang'],p_dict['maxout'],nextra,ntotal,qextra,qtotal,p_dict['min_amp'],len(f_p_pol))
        faultnorms_all=faultnorms_all[:,:nf]
        faultslips_all=faultslips_all[:,:nf]
        strike_all=strike_all[:nf]
        dip_all=dip_all[:nf]
        rake_all=rake_all[:nf]
    else: # Python version
        faultnorms_all,faultslips_all=fun.focal_gridsearch(sr_azimuth,takeoff,p_pol,sp_amp,dir_cos_dict,nextra,ntotal,qextra,qtotal,p_dict['maxout'],dir_cos_dict['ncoor'])

        # Calculates strike,dip,rake from normal,slip vectors for output
        if ((len(p_dict['outfile2'])>0) | (p_dict['plot_acceptable_solutions'])):
            strike_all,dip_all,rake_all=fun.sdr_from_vector(faultnorms_all,faultslips_all)

    if (faultnorms_all).shape[1]==0:
        print('{} / {}\t({})\n'.format(event_x,num_events-1,event_id)+
              '\tNo solution found for {}'.format(event_id))
        # return 0
        return mech_dict

    p_the_mc=takeoff[:,0]
    p_azi_mc=sr_azimuth[:,0]

    # Calculates the probabilities for potential mech solutions
    mech_df=fun.mech_probability(faultnorms_all,faultslips_all,p_dict['cangle'],p_dict['prob_max'],iterative_avg=p_dict['iterative_avg'])

    if len(mech_df)==0: # No accepted solution
        print('{} / {}\t({})\n'.format(event_x,num_events-1,event_id)+
              '\tNo accepted solution found for {}'.format(event_id))
        # return 0
        return mech_dict

    # Calculates the misfit for prefered solutions
    mech_df,pol_agreement_out,sp_diff_out=fun.mech_misfit(mech_df,p_azi_mc,p_the_mc,p_pol,sp_amp)

    # Mech solution quality rating
    mech_df=fun.mech_quality(mech_df,qual_criteria_dict)
    mech_df['num_p_pol']=np.sum(p_pol!=0)
    mech_df['num_sp_ratios']=np.sum(~np.isnan(sp_amp))

    if p_dict['min_quality_report']:
        mech_df=mech_df.loc[mech_df['qual']<=p_dict['min_quality_report'],:].reset_index(drop=True)
        if len(mech_df)==0: # No accepted solution
            print('{} / {}\t({})\n'.format(event_x,num_events-1,event_id)+
                '\tNo solution met the minimum quality ({}) for {}'.format(p_dict['min_quality_report'],event_id))
            # continue
            # return 0
            return mech_dict

    # Rounds mech solution values
    angle_col=['str_avg', 'dip_avg', 'rak_avg','rms_diff','rms_diff_aux']
    quality_col=['prob','mfrac','mavg','stdr']
    mech_df[angle_col]=mech_df[angle_col].round(p_dict['output_angle_precision'])
    mech_df[quality_col]=(mech_df[quality_col]*100).round(p_dict['output_quality_precision'])

    # Writes preferred mechanisms to file
    if p_dict['outfile1']:
        if len(cat_df):
            event_df=cat_df.loc[cat_df.event_id==event_id].iloc[0]
        else:
            event_df=[]
        out.write_outfile1(p_dict['outfile1'],mech_df,event_df,event_id)

    # Writes acceptable mechanisms to file
    if p_dict['outfile2']:
        out.write_outfile2(p_dict['outfile2'],event_id,strike_all,dip_all,rake_all,faultnorms_all,faultslips_all,p_dict['output_angle_precision'],p_dict['output_vector_precision'])

    if (p_dict['stfile']) and (p_dict['outfile_pol_info']):
        out_takeoff=np.round(p_the_mc,p_dict['output_angle_precision'])
        out_sr_az=np.round(p_azi_mc,p_dict['output_angle_precision'])
    else:
        out_takeoff=-999
        out_sr_az=-999

    # Plots the focal mechanism and saves the figure
    if p_dict['outfolder_plots']:
        import functions.plot_mech as plot_mech # For plotting mechanism solutions
        if p_dict['plot_acceptable_solutions']:
            plot_mech.plot_mech(mech_df,event_pol_df,takeoff[:,0],sr_azimuth[:,0],p_dict,acceptable_sdr=np.vstack([strike_all,dip_all,rake_all]))
        else:
            plot_mech.plot_mech(mech_df,event_pol_df,takeoff[:,0],sr_azimuth[:,0],p_dict)

    # Calculates the stdev of the takeoff and azimuths
    if p_dict['outfile_pol_info'] and p_dict['stfile']:
        takeoff_uncertainty_out=np.std(takeoff,axis=1).round(p_dict['output_angle_precision'])

        rad_sr_azimuth=np.deg2rad(sr_azimuth)
        sr_sin_mean=np.mean(np.sin(rad_sr_azimuth),axis=1)
        sr_cos_mean=np.mean(np.cos(rad_sr_azimuth),axis=1)
        R=np.minimum(1,np.hypot(sr_sin_mean, sr_cos_mean))
        azimuth_uncertainty_out=np.degrees(np.sqrt(-2*np.log(R))).round(p_dict['output_angle_precision'])
    else:
        takeoff_uncertainty_out=-1.
        azimuth_uncertainty_out=-1.

    # Prints summary results for the event to std out
    print('{} / {}\t({})\n\tS: {}   D: {}   R: {}   U: {}   Q: {}\n\tRuntime: {:.2f} sec'.format(
                    event_x,
                    num_events-1,
                    event_id,
                    mech_df.loc[0,'str_avg'],
                    mech_df.loc[0,'dip_avg'],
                    mech_df.loc[0,'rak_avg'],
                    mech_df.loc[0,'rms_diff'],
                    mech_df.loc[0,'qual'],
                    (time.time()-event_runtime_start) ),flush=True)

    # return event_pol_df
    return {'event_index':event_pol_df.index,'pol_agreement_out':pol_agreement_out,
            'sp_diff_out':sp_diff_out,'takeoff':out_takeoff,'sr_az':out_sr_az,
            'takeoff_uncertainty':takeoff_uncertainty_out,'azimuth_uncertainty':azimuth_uncertainty_out,
            'mech_qual':mech_df.loc[0,'qual']}
