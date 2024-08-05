import os
import pandas as pd

# Define a class to run SKHASH with new polarities
class RerunSKHASH:
    def __init__(self, **kwargs):
        print(kwargs.keys())
        # Unpack common kwargs
        self.pick_pol_path = kwargs.get('pick_pol_path', None)
        self.conpfile     = kwargs.get('pol_path', None)
        self.skhash_dir         = kwargs.get('skhash_dir', 'SKHASH2/SKHASH')

        # Generate SKHASH format polarity file
        pol_file = self.csv_picks2skhash_pol_file(self.pick_pol_path
                    ).to_csv(self.conpfile, index=False)

        # Edit the SKHASH control file
        self.edit_skhash_control_file(**kwargs)

        # Run SKHASH
        self.make_run_skhash_python_script(**kwargs)


    def csv_picks2skhash_pol_file(self, picks_path, outpath=None, write=False):
        """
        Convert the NCEDC-style picks csv file to SKHASH format.

        Input: 
            NCEDC-style picks csv file -
            event_id,etime,elat,elon,edep,emag,station_id,phase_type,phase_polarity,phase_time,sta,slon,slat,sta-event_dist_km,velocity_kmps,SNR

        Output: SKHASH polarity format csv file -
            event_id, event_id2, station, location, channel, p_polarity, origin_latitude, origin_longitude, origin_depth_km
        """

        # Read picks and filter for P picks and non-missing polarities
        picks_df = pd.read_csv(picks_path)
        picks_df = picks_df[picks_df['phase_type'] == 'P'].dropna(subset=['phase_polarity'])

        # Define a dictionary for polarity conversion
        polarity_map = {'U': 1, 'D': -1, 'unknown': 0}

        # Create SKHASH format dataframe
        df = picks_df[['event_id']].copy()
        df['event_id2'] = picks_df['event_id']
        df['station'] = picks_df.station_id.astype(str).str.split('.').str[1]
        df['location'] = picks_df['station_id'].str.split('.', expand=True)[2].apply(lambda x: '--' if x == '' else x)
        df['channel'] = picks_df['station_id'].str.split('.', expand=True)[3].apply(lambda x: x+'Z' if len(x) == 2 else x)
        df['p_polarity'] = picks_df['phase_polarity'].replace(polarity_map)
        df[['origin_latitude', 'origin_longitude', 'origin_depth_km']] = picks_df[['elat', 'elon', 'edep']].astype(float).round(4)

        if write:
            df.to_csv(outpath, index=False)
        else:
            return df


    def make_run_skhash_python_script(self, **kwargs):

        control_file_path = kwargs.get('control_file_path')
        skhash_dir        = kwargs.get('skhash_dir')

        # Change the current working directory to skhash_dir
        os.chdir(skhash_dir)

        # Run SKHASH
        python_script = f"python SKHASH.py {control_file_path}"

        # Run the script
        print("Running the script...")
        # os.system(f'conda activate obspy | {python_script}')
        os.system(python_script)

        
    def edit_skhash_control_file(self, **kwargs):
        """
        Edit the SKHASH control file to include the paths to the station, polarity, and reverse files.
        """

        # Unpack the kwargs
        output_path = kwargs.get('control_file_path', None)
        # input files
        conpfile = kwargs.get('pol_path', None)
        stfile = kwargs.get('stn_path', None)
        # output files
        outfile1 = kwargs.get('mech_path', None)
        outfile2 = kwargs.get('alt_mech_path', None)
        outfile_pol_agree = kwargs.get('pol_agree_path', None)
        outfile_pol_info = kwargs.get('pol_info_path', None)
        # control parameters
        vmodel_paths = kwargs.get('vmodel_paths', None)
        output_angle_precision = kwargs.get('output_angle_precision', 4)
        require_network_match = kwargs.get('require_network_match', False)
        allow_duplicate_stations = kwargs.get('allow_duplicate_stations', True)
        min_polarity_weight = kwargs.get('min_polarity_weight', 0)
        dang = kwargs.get('dang', 5)
        nmc = kwargs.get('nmc', 50)
        maxout = kwargs.get('maxout', 500)
        ratmin = kwargs.get('ratmin', 2)
        badfrac = kwargs.get('badfrac', 0.1)
        qbadfrac = kwargs.get('qbadfrac', 0.3)
        delmax = kwargs.get('delmax', 175)
        max_pgap = kwargs.get('max_pgap', 65)
        cangle = kwargs.get('cangle', 45)
        prob_max = kwargs.get('prob_max', 0.2)
        azmax = kwargs.get('azmax', 0)
        max_agap = kwargs.get('max_agap', 135)
        outfolder_plots = kwargs.get('outfolder_plots', None)
        plot_station_names = kwargs.get('plot_station_names', False)
        plot_acceptable_solutions = kwargs.get('plot_acceptable_solutions', False)

        # Create the control file
        control_file = f"""## Control file for SKHASH
$input_format  # format of input files
SKHASH

$num_cpus      # number of CPUs to use
0

$conpfile        # P-polarity input filepath
{conpfile}

$stfile        # station list filepath
{stfile}

$outfile1      # focal mechanisms output filepath 
{outfile1}

$outfile2 : # Path to acceptable plane output file
{outfile2}

$outfile_pol_agree  # record of polarity (dis)agreeement output filepath # examples/maacama_SKHASH_MTJ/OUT/out_polagree.csv
{outfile_pol_agree}

$outfile_pol_info # examples/maacama_SKHASH_MTJ/OUT/out_polinfo.csv
{outfile_pol_info}

$vmodel_paths  # whitespace/newline delimited list of paths to the velocity models 
{vmodel_paths}

$output_angle_precision
{output_angle_precision}

$require_network_match
{require_network_match}

$allow_duplicate_stations
{allow_duplicate_stations}

$min_polarity_weight
{min_polarity_weight}

$dang          # minimum grid spacing (degrees)
{dang}

$nmc           # number of trials (e.g., 30)
{nmc}

$maxout        # max num of acceptable focal mech. outputs (e.g., 500)
{maxout}

$ratmin        # minimum allowed signal to noise ratio
{ratmin}

$badfrac       # fraction polarities assumed bad
{badfrac}

$qbadfrac      # assumed noise in amplitude ratios, log10 (e.g. 0.3 for a factor of 2)
{qbadfrac}

$delmax        # maximum allowed source-receiver distance in km.
{delmax}

$max_pgap      # maximum allowed takeoff angle in degrees
{max_pgap}

$cangle        # angle for computing mechanisms probability
{cangle}

$prob_max      # probability threshold for multiples (e.g., 0.1)
{prob_max}

$azmax         # Maximum allowed source-station azimuth uncertainty in degrees [0 = all allowed]
{azmax}

$max_agap      # maximum azimuthal gap in degrees
{max_agap}

$allow_hypocenters_outside_table # False: only hypocenters within the velocity model are allowed
True

"""
    
        if outfolder_plots is not None:
            control_file += f"""\n$outfolder_plots : #Path to folder where simple focal mechanism plots \n{outfolder_plots}\n"""
        if plot_station_names:
            control_file += f"""\n$plot_station_names : #Plot station names \n{plot_station_names}\n"""
        if plot_acceptable_solutions:
            control_file += f"""\n$plot_acceptable_solutions : #Plot acceptable solutions \n{plot_acceptable_solutions}\n"""

        # Write the control file
        if output_path is not None:
            with open(output_path, 'w') as f:
                print(f"Writing the control file to\n {output_path}")
                f.write(control_file)
        else:
            return control_file


# Example usage
if __name__ == "__main__":
    RerunSKHASH()