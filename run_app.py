import os
from PyQt5.QtWidgets import QApplication

# Import sef-defined modules
from fun_beachball_widget import MainApp

user = 'mdarifulislam'

# Define the paths (Use absolute paths in most cases)
file_dir      = f'/Users/{user}/Library/CloudStorage/OneDrive-IndianaUniversity/Research/Github/FM2STRESS/FM2STRESS_project/data/NCEDC_picks/HASH_IN_OUT'
mseed_dir     = f'/Users/{user}/Library/CloudStorage/OneDrive-IndianaUniversity/Research/Data/NCEDC_events_data/mseed/2008'
mech_path     = os.path.join(file_dir, 'OUT/out_NC_PN_QC_2008.csv')
alt_mech_path = os.path.join(file_dir, 'OUT/out2_NC_PN_QC_2008.csv')
pol_info_path = os.path.join(file_dir, 'OUT/out_polinfo_NC_PN_QC_2008.csv')
pol_agree_path= os.path.join(file_dir, 'OUT/out_polagree_NC_PN_QC_2008.csv') 
pol_path      = os.path.join(file_dir, 'IN/pol_concensus_NC_PN_QC_2008.csv')
stn_path      = os.path.join(file_dir, 'IN/station_master.csv')
pick_pol_path = f'/Users/{user}/Library/CloudStorage/OneDrive-IndianaUniversity/Research/Data/NCEDC_events_data/Markers/2008_pyrocko_markers_filt.csv'
skhash_root_dir = os.path.abspath('FM2STRESS_project/code/InteractiveFM/SKHASH2/SKHASH')

# Make app and window for beachball plot
app = QApplication([])

input_params = {
    'mech_path': mech_path,         # (outfile1) SKHASH results
    'alt_mech_path': alt_mech_path, # (outfile2) SKHASH output accepted mechanisms
    'pol_path': pol_path,           # SKHASH input polarities
    'pol_info_path': pol_info_path, # (outfile_pol_info) SKHASH output polarities (same as pol_path with additional info)
    'pol_agree_path': pol_agree_path,# (outfile_pol_agree) SKHASH output polarities with agreement
    'stn_path': stn_path,           # Master station csv (SKHASH format)
    'pick_pol_path': pick_pol_path, # phase picks and polarities csv (converted from pyrocko markers, filtered)
    'mseed_dir': mseed_dir,
    'hor_line': True,
    'zoom': 2,
    'slice_len': 0.5,
    'normalize': True,
    # SKHASH control file parameters
    'vmodel_paths': os.path.join(skhash_root_dir, 'examples/velocity_models_MTJ/vz_MTJ.txt'),
    'max_agap': 170,
    'delmax': 0,
    # Rerun SKHASH parameters
    'mini_or_ana': 'miniconda3',
    'skhash_dir': skhash_root_dir,
    'control_file_path': os.path.join(file_dir, 'control_file_app.txt')
}

# Make app and window for beachball and waveform plot
main_window = MainApp(**input_params)

# Show the main window
main_window.show()

# Execute the application
app.exec_()