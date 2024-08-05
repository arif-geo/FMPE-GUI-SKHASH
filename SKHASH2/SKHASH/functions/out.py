# External libraries
import numpy as np
import pandas as pd


def create_outfile1(outfile1,cat_df,pol_df):
    '''
    Creates the outfile of the preferred mech solutions
    '''
    header_str_out='event_id,strike,dip,rake,quality,fault_plane_uncertainty,aux_plane_uncertainty,num_p_pol,num_sp_ratios,'+\
                        'polarity_misfit,prob_mech,sta_distribution_ratio,sp_misfit,mult_solution_flag'
    if len(cat_df):
        header_str_out+=','
        if 'origin_DateTime' in cat_df:
            header_str_out+='time,'
        if 'event_mag' in pol_df:
            header_str_out+='magnitude,'
        header_str_out+='origin_lat,origin_lon,origin_depth_km,horz_uncert_km,vert_uncert_km'
    with open(outfile1, "w") as f_outfile1:
        f_outfile1.write(header_str_out+'\n')
    return True


def create_outfile2(outfile2):
    '''
    Creates the outfile of the acceptable mech solutions
    '''
    with open(outfile2, "w") as f_outfile2:
        f_outfile2.write('event_id,mech_number,strike,dip,rake,norm_N,norm_E,norm_Z,norm_aux_N,norm_aux_E,norm_aux_Z\n')
    return True


def write_outfile1(outfile1,mech_df,event_df,event_id):
    '''
    Writes mech solutions for an event to outfile1.
    '''
    event_str_output=''
    if len(event_df):
        event_str_output+=','
        if 'origin_DateTime' in event_df:
            event_str_output+=event_df.origin_DateTime.strftime('%Y-%m-%d %X')+','
        if 'event_mag' in event_df:
            event_str_output+=str(event_df['event_mag'])+','
        event_str_output+=str(event_df['origin_lat'])+','
        event_str_output+=str(event_df['origin_lon'])+','
        event_str_output+=str(event_df['origin_depth_km'])+','
        event_str_output+=str(event_df['horz_uncert_km'])+','
        event_str_output+=str(event_df['vert_uncert_km'])

    with open(outfile1, "a") as f_outfile1:
        if len(mech_df)>1:
            mech_df['mflag']=True
        else:
            mech_df['mflag']=False
        for imult in range(len(mech_df)):
            f_outfile1.write(('{},{},{},{},{},'+
                              '{},{},{},{},{},'+
                              '{},{},{},{}{}\n').format(
                event_id,
                mech_df.loc[imult,'str_avg'], # strike
                mech_df.loc[imult,'dip_avg'], # dip
                mech_df.loc[imult,'rak_avg'], # rake
                mech_df.loc[imult,'qual'], # ad-hoc mech quality
                mech_df.loc[imult,'rms_diff'], # fault plane uncertainty
                mech_df.loc[imult,'rms_diff_aux'], # aux plane uncertainty
                mech_df.loc[imult,'num_p_pol'], # num p polarity picks
                mech_df.loc[imult,'num_sp_ratios'], # num S/P ratios
                mech_df.loc[imult,'mfrac'], # weighted percent misfit of first motions
                mech_df.loc[imult,'prob'], # probability mechanism close to solution
                mech_df.loc[imult,'stdr'], # 100*(station distribtuion ratio)
                mech_df.loc[imult,'mavg'], # 100*(average log10(S/P) misfit)
                mech_df.loc[imult,'mflag'], # Flag indicating whether there are multiple solutions for the event
                event_str_output
                ))
    return True


def write_outfile2(outfile2,event_id,strike_all,dip_all,rake_all,faultnorms_all,faultslips_all,output_angle_precision,output_vector_precision):
    '''
    Writes acceptable mech solutions for an event to outfile2.
    '''
    acceptable_mechs_df=pd.DataFrame(columns=['event_id','mech_number','strike','dip','rake',
                                              'norm_N','norm_E','norm_Z','norm_aux_N','norm_aux_E','norm_aux_Z'],
                                              index=np.arange(len(strike_all)))

    acceptable_mechs_df['event_id']=event_id
    acceptable_mechs_df['mech_number']=np.arange(len(acceptable_mechs_df))
    acceptable_mechs_df['strike']=np.round(strike_all,output_angle_precision)
    acceptable_mechs_df['dip']=np.round(dip_all,output_angle_precision)
    acceptable_mechs_df['rake']=np.round(rake_all,output_angle_precision)
    acceptable_mechs_df['norm_N']=np.round(faultnorms_all[0,:],output_vector_precision)
    acceptable_mechs_df['norm_E']=np.round(faultnorms_all[1,:],output_vector_precision)
    acceptable_mechs_df['norm_Z']=np.round(faultnorms_all[2,:],output_vector_precision)
    acceptable_mechs_df['norm_aux_N']=np.round(faultslips_all[0,:],output_vector_precision)
    acceptable_mechs_df['norm_aux_E']=np.round(faultslips_all[1,:],output_vector_precision)
    acceptable_mechs_df['norm_aux_Z']=np.round(faultslips_all[2,:],output_vector_precision)
    acceptable_mechs_df.to_csv(outfile2,mode='a',sep=',',index=False,header=False)
    return True


def pol_agree(pol_df,p_dict):
    '''
    Writes the polarity agreement records for all events to outfile_pol_agree
    '''
    if (p_dict['min_quality_report']):
        pol_df=pol_df[pol_df.mech_quality<=p_dict['min_quality_report']]

    if not(pol_df.empty):
        pol_df['pol_agreement_weighted']=pol_df['pol_agreement']*(pol_df['p_polarity'].abs())

        group_pol_df=pol_df.groupby(['sta_code','source'])

        pol_agree_df=group_pol_df.agg({'pol_agreement':['sum','count'],'p_polarity':[lambda lam: (abs(lam).sum())],'pol_agreement_weighted':['sum'] })
        pol_agree_df['count_correct']=pol_agree_df['pol_agreement']['sum'].astype(int)
        pol_agree_df['count_total']=pol_agree_df['pol_agreement']['count'].astype(int)
        pol_agree_df['pol_accuracy']=(100*pol_agree_df['count_correct']/pol_agree_df['count_total']).round(1)

        pol_agree_df['weight_correct']=pol_agree_df['pol_agreement_weighted']['sum'].astype(float)
        pol_agree_df['weight_total']=pol_agree_df['p_polarity']['<lambda>'].astype(float)
        pol_agree_df['weight_pol_accuracy']=(100*pol_agree_df['weight_correct']/pol_agree_df['weight_total']).round(1)
        pol_agree_df['weight_correct']=pol_agree_df['weight_correct'].round(5)
        pol_agree_df['weight_total']=pol_agree_df['weight_total'].round(5)

        pol_agree_df=pol_agree_df.drop(columns=['pol_agreement','p_polarity','pol_agreement_weighted'])
        pol_agree_df=pol_agree_df.reset_index()
        pol_agree_df.columns=pol_agree_df.columns.droplevel(1)
        pol_agree_df=pol_agree_df.sort_values(by=['weight_pol_accuracy','count_correct','sta_code'],ascending=[True,False,True])

        pol_agree_df=pol_agree_df.loc[pol_agree_df.count_total>0,:].reset_index(drop=True)

        pol_agree_df.to_csv(p_dict['outfile_pol_agree'],index=False)
    return True


def sp_agree(pol_df,p_dict):
    '''
    Writes the S/P agreement records for all events to outfile_sp_agree
    '''
    if (p_dict['min_quality_report']):
        pol_df=pol_df[pol_df.mech_quality<=p_dict['min_quality_report']]
    if not(pol_df.empty):

        group_sp_agree_df=pol_df.groupby(['sta_code','source'])

        sp_agree_df=group_sp_agree_df['sp_diff'].agg(['mean','std','count']).reset_index()
        sp_agree_df=sp_agree_df.drop(np.where(sp_agree_df['count']==0)[0])

        sp_agree_df=sp_agree_df.sort_values(by=['mean','count','std'],ascending=[True,False,True])
        sp_agree_df.loc[pd.isnull(sp_agree_df['std']),'std']=-999.

        sp_agree_df['mean']=(sp_agree_df['mean']).round(4)
        sp_agree_df['std']=(sp_agree_df['std']).round(4)
        sp_agree_df=sp_agree_df.rename(columns={'mean':'avg_sp_diff','std':'std_sp_diff'})

        sp_agree_df.to_csv(p_dict['outfile_sp_agree'],index=False)
    return True


def pol_info(pol_df,p_dict):
    '''
    Writes information about the pol info for all events to outfile_pol_info
    '''
    if 'sp_ratio' in pol_df.columns:
        pol_df['sp_ratio']=(10**pol_df['sp_ratio']).round(5)
    if 'sp_diff' in pol_df.columns:
        pol_df['sp_diff']=(10**pol_df['sp_diff']).round(5)

    if not(pol_df.empty):
        pol_df.to_csv(p_dict['outfile_pol_info'],index=False)
    return True
