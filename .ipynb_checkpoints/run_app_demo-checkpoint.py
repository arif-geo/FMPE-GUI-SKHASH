import os
from PyQt5.QtWidgets import QApplication

# Import sef-defined modules
from fun_beachball_widget import MainApp

# Define the paths (Use absolute paths in most cases)
file_dir      = os.path.abspath('./TestData')
# output files (will be created if not provided)
mech_path     = os.path.join(file_dir, 'skhash_output/out_mech.csv')
alt_mech_path = os.path.join(file_dir, 'skhash_output/out2_alt_mech.csv')
pol_info_path = os.path.join(file_dir, 'skhash_output/out_polinfo.csv')
pol_agree_path= os.path.join(file_dir, 'skhash_output/out_polagree.csv') 
# input files
pol_path      = os.path.join(file_dir, 'skhash_input/pol_concensus.csv')   # will be created if not provided
stn_path      = os.path.join(file_dir, 'skhash_input/Station_master_skhash_format.csv')
pick_pol_path = os.path.join(file_dir, 'skhash_input/PhasePicks_my_standard.csv')
eq_cat_path   = os.path.join(file_dir, 'skhash_input/NCEDC_eq_cat_above_slab.csv')
mseed_dir     = os.path.join(file_dir, 'mseed')
skhash_root_dir = os.path.abspath('SKHASH2/SKHASH')

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
    'eq_cat_path': eq_cat_path,     # Earthquake catalog csv
    'mseed_dir': mseed_dir,
    'hor_line': True,
    'zoom': 2,
    'slice_len': 0.5,
    'normalize': True,
    # beachball plot parameters
    'acceptable_sdr': True,
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
