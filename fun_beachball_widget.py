import os
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtCore import Qt, pyqtSignal
# from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QVBoxLayout, QWidget, QSplitter, QFrame
from PyQt5.QtWidgets import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

from obspy import read, UTCDateTime

# Custom modules
from fun_plot_beachball import plot_mech, get_nearest_station, get_exact_station
from fun_makefile_run_SKHASH import RerunSKHASH

class BeachballPlot(QWidget):
    clicked_station_signal = pyqtSignal(str, str)  # Define signal for clicked station

    def __init__(self, mech_df=None, alt_mech_df=None, pol_df=None, pol_info_df=None, **kwargs):
        super().__init__()

        self.mech_df = mech_df
        self.alt_mech_df = alt_mech_df
        self.pol_df = pol_df
        self.pol_info_df = pol_info_df

        # Plot the mechanism using pre-defined data
        fig, ax, sc_up, sc_down, sc_zero, xy_data, event_id = plot_mech(
            self.mech_df, self.pol_df, self.pol_info_df, self.alt_mech_df, **kwargs)

        self.fig = fig
        self.ax = ax
        self.sc_up = sc_up
        self.sc_down = sc_down
        self.sc_zero = sc_zero
        self.xy_data = xy_data
        self.event_id = event_id

        # Connect pick event to onpick function with signal emission
        self.fig.canvas.mpl_connect('pick_event', self.onpick)

        layout = QVBoxLayout(self)
        self.canvas = FigureCanvasQTAgg(self.fig)
        layout.addWidget(self.fig.canvas)

        # Add the navigation toolbar
        self.toolbar = NavigationToolbar2QT(self.canvas, self)
        layout.addWidget(self.toolbar)
        self.setLayout(layout)

    def onpick(self, event):
        '''
        Station_name full station name (e.g. 'NC.PFO..EHZ' not just 'PFO')
        '''
        if event.artist in [self.sc_up, self.sc_down, self.sc_zero]:  # Check if clicked on beachball collection
            point = event.artist.get_offsets()[event.ind[0]]
            # clicked_station = get_exact_station(point, self.xy_data, self.pol_info_df)
            clicked_station = get_nearest_station(point, self.xy_data, self.pol_info_df)
            self.clicked_station_signal.emit(self.event_id, clicked_station)

# Lets define another class which will replace the SimpleApp class and plot seismic waveforms
class WaveformPlot(QWidget):
    def __init__(self):
        super().__init__()
        self.fig, self.ax = plt.subplots()
        self.canvas = self.fig.canvas
        
        layout = QVBoxLayout(self)
        layout.addWidget(self.canvas)
        self.setLayout(layout)

    def plot_waveform(self, event_id, station_name, pick_pol_df, data_dir, **kwargs):

        # Optional parameters
        slice_len = kwargs.get('slice_len', 0.5)    # if not provided, default is 0.5
        zoom = kwargs.get('zoom', 1)                # if not provided, default is 1
        normalize = kwargs.get('normalize', True)   # if not provided, default is True
        hor_line = kwargs.get('hor_line', True)     # if not provided, default is True

        if event_id and station_name and data_dir and pick_pol_df is not None:
            # Phase pick data
            pick_df = pick_pol_df.loc[(pick_pol_df['event_id'] == event_id) & 
                                        (pick_pol_df['station_id'].str.split('.').str[1] == station_name.split('.')[0])]
            if len(pick_df) > 1:
                pick_df = pick_df.iloc[0]
            pick_df_idx = pick_df.index.values[0]
            phase_time = UTCDateTime(pd.to_datetime(pick_df['phase_time'].values[0]))

            # Load waveform data
            st = read(f"{data_dir}/{event_id}.mseed")
            tr = st.select(station=station_name.split('.')[0])[0]
            tr.detrend('demean')
            tr.taper(0.001)
            tr.filter('bandpass', freqmin=1.5, freqmax=10)
            tr = tr.slice(starttime=phase_time-slice_len, endtime=phase_time+slice_len)
            if normalize: 
                tr.normalize()
            if zoom:
                tr.data = tr.data * zoom
            xtime = np.arange(0, tr.stats.npts / tr.stats.sampling_rate, tr.stats.delta)

            self.ax.clear()
            self.ax.plot(xtime, tr.data)                                # Plot waveform
            self.ax.axvline(x=phase_time - tr.stats.starttime, color='r', linestyle='--', lw=0.5)    # plot vertical line at phase pick time
            self.ax.text(phase_time - tr.stats.starttime - 0.5, 0.5, f"{tr.id}\nPol:{pick_df['phase_polarity'].values[0]}\nIdx{pick_df_idx}", fontsize=8)
            if hor_line:
                self.ax.axhline(y=0, color='k', linestyle='--', lw=0.5)
            self.ax.set_title(f'Waveform for Event {event_id}, Station {station_name}')
            self.ax.set_ylim(-1.5, 1.5)
            self.ax.set_xlabel('Time (s)')
            self.ax.set_ylabel('Amplitude [Normalized]')
            # tight_layout() 
            self.fig.tight_layout()
        else:
            # Empty plot
            self.ax.clear()
            self.ax.set_title('No Data')
        self.canvas.draw()
    
    def update_waveform(self, event_id, station_name, pick_pol_df, data_dir, **kwargs):
        self.plot_waveform(event_id, station_name, pick_pol_df, data_dir, **kwargs)
        


class MainApp(QMainWindow):
    '''

    Description of the input parameters:
    - mech_path: CSV containing focal mechanism data (SKHASH results)
    - alt_mech_path: CSV containing accepted focal mechanisms (SKHASH output)
    - pol_path: CSV  containing polarity data (SKHASH input)
    - pol_info_path: CSV containing polarity data with additional information (SKHASH output)
    - stn_path: CSV containing station information (SKHASH format)
    - *** pick_pol_path: (Very crucial, will be used to rerun SKHASH with new polarities)
        CSV containing phase picks and polarities
        (Converted from pyrocko markers following PhaseNet format,
        preferably filter the data in a separate script before passing here,
        fillers include:
            - anomolous high velocity picks,
            - Magnitude based distance filter [e.g. M<1.8 with 80km radius],
            - SNR based filter [e.g. SNR>10 at least >6 for good picks],
            - Column names follow PhaseNet format but not exactly the same:
            "event_id", "etime", "elat", "elon","edep","emag","station_id","phase_type",
            "phase_polarity","phase_time","sta","slon","slat","sta-event_dist_km","velocity_kmps"
            and so on...)

    '''
    def __init__(self, **kwargs):
        super().__init__()
        
        # Unpack the keyword arguments
        self.mech_path = kwargs.get('mech_path', None)
        self.alt_mech_path = kwargs.get('alt_mech_path', None)
        self.pol_path = kwargs.get('pol_path', None)
        self.pol_info_path = kwargs.get('pol_info_path', None)
        self.stn_path = kwargs.get('stn_path', None)
        self.pick_pol_path = kwargs.get('pick_pol_path', None)
        self.mseed_dir = kwargs.get('mseed_dir', None)

        # Initialize event_id and station_name to None
        self.event_id = None
        self.station_name = None
        self.current_event_index = 0
        
        # Load all data to start with
        self.initUI(**kwargs)
        self.load_data(**kwargs)
    

    def initUI(self, **kwargs):
        # Create main window and layout
        self.setWindowTitle('Beachball and Waveform Plot')
        self.setGeometry(0, 0, 1500, 800)
        self.central_widget = QWidget()
        self.layout = QVBoxLayout(self.central_widget)
        self.setCentralWidget(self.central_widget)

        # Create splitter and frames
        self.splitter = QSplitter(Qt.Horizontal)
        self.beachball_frame = QFrame()
        self.waveform_frame = QFrame()
        self.button_frame = QFrame()

        # Create button layout
        self.button_layout = QVBoxLayout()
        self.button_frame.setLayout(self.button_layout)

        # Add buttons: Reverse polarity
        self.add_button_with_label(
            'Reverse polarity \n(i.e. -1 will be 1 and vice versa)', 
            'Reverse Polarity', self.reverse_polarity)

        # Add buttons: Set polarity to zero
        self.add_button_with_label(
            'Set Polarity to Zero for this station', 'Set Polarity to Zero', self.set_polarity_zero)

        # Add buttons: Next Event
        self.add_button_with_label('Go to Next/Previous Event', 'Next Event >>', self.load_next_event, **kwargs)
        self.add_button_with_label('', '<< Previous Event', self.load_next_event, reverse=True, **kwargs)

        # Add buttons: Save edited polarity file
        self.add_button_with_label('Save Edited Polarity File', 'Save Polarity File', self.save_polarity_file)

        # Add buttons: Rerun SKHASH with new polarities
        self.add_button_with_label(
            'Re-run SKHASH with new polarities\n [Must save the edited polarity file first]', 
            'Re-run SKHASH', self.rerun_skhash, **kwargs)

        # Add splitter to main layout
        self.layout.addWidget(self.splitter)

        # Add frames to splitter
        self.splitter.addWidget(self.beachball_frame)
        self.splitter.addWidget(self.button_frame)
        self.splitter.addWidget(self.waveform_frame)

    def update_ui(self, **kwargs):

        # Create new layouts for the frames
        beachball_layout = QVBoxLayout()
        waveform_layout = QVBoxLayout()

        # Create the beachball and waveform plots
        self.bb_plot = self.plot_beachball_widget(self.mech_df, self.alt_mech_df, self.pol_df, self.pol_info_df, **kwargs)
        self.wf_plot = WaveformPlot() # Inputs for this plot will be added in Func 'plot_beachball_widget' by calling 'store_event_station'

        # Add plots to their respective layouts
        beachball_layout.addWidget(self.bb_plot)
        waveform_layout.addWidget(self.wf_plot)

        # Remove existing layouts if they exist
        if self.beachball_frame.layout() is not None:
            QWidget().setLayout(self.beachball_frame.layout())
        if self.waveform_frame.layout() is not None:
            QWidget().setLayout(self.waveform_frame.layout())

        # Set the layouts for the frames
        self.beachball_frame.setLayout(beachball_layout)
        self.waveform_frame.setLayout(waveform_layout)
            
        # Update the splitter to reflect changes
        self.splitter.update()


    def add_button_with_label(self, label_text, button_text, callback, **kwargs):
        '''
        - label_text: text to be displayed next to the button
        - button_text: text to be displayed on the button
        - callback: function to be called when the button is clicked
        '''
        layout = QVBoxLayout()
        label = QLabel(label_text)
        button = QPushButton(button_text, self)
        button.clicked.connect(lambda: callback(**kwargs))
        if not label_text == '':
            layout.addWidget(label, alignment=Qt.AlignCenter)
        layout.addWidget(button)
        self.button_layout.addLayout(layout)

    def load_data(self, **kwargs):
        self.stn_df_all = pd.read_csv(self.stn_path)
        self.pick_pol_df_all = pd.read_csv(self.pick_pol_path)
        # If mech_df data does not exist
        if os.path.exists(self.mech_path):
            self.mech_df_all = pd.read_csv(self.mech_path)
            self.alt_mech_df_all = pd.read_csv(self.alt_mech_path)
            self.pol_info_df_all = pd.read_csv(self.pol_info_path)
            self.pol_df_all = pd.read_csv(self.pol_path)

            if kwargs.get('eq_cat_path'): # If earthquake catalog is provided
                self.eq_df_all = pd.read_csv(kwargs.get('eq_cat_path', None), parse_dates=['time'])
                # Merge by event_id to add time column
                self.mech_df_all = pd.merge(
                    self.mech_df_all, self.eq_df_all[['id', 'time']], left_on='event_id', right_on='id', how='left'
                    ).drop(columns='id').sort_values(by='time')

            # Get the unique event IDs
            self.event_ids = self.mech_df_all['event_id'].unique()
            self.load_event_data(self.current_event_index, **kwargs)

        else:
            # Only run SKHASH to get the mechanism data
            self.initUI(**kwargs)
            

    def load_event_data(self, event_index, **kwargs):
        self.event_id = self.event_ids[event_index]

        self.mech_df = self.mech_df_all.loc[self.mech_df_all['event_id'] == self.event_id].iloc[0:1]
        self.alt_mech_df = self.alt_mech_df_all.loc[self.alt_mech_df_all['event_id'] == self.event_id]
        self.pol_info_df = self.pol_info_df_all.loc[self.pol_info_df_all['event_id'] == self.event_id]
        self.pol_df = self.pol_df_all.loc[self.pol_df_all['event_id'] == self.event_id]
        self.pick_pol_df = self.pick_pol_df_all.loc[self.pick_pol_df_all['event_id'] == self.event_id]
        
        # Update the UI with the new event data  
        self.update_ui(**kwargs)
        
    def load_next_event(self,reverse=False, **kwargs):
        # clear matplotlib figure
        plt.close('all')
        
        # Update the current event index
        if reverse: # Previous event
            self.current_event_index = (self.current_event_index - 1) % len(self.event_ids)
        else: # Next event
            self.current_event_index = (self.current_event_index + 1) % len(self.event_ids)

        # Load the data for the next event
        self.load_event_data(self.current_event_index, **kwargs)

    # define a separate func to plot beachball
    def plot_beachball_widget(self, mech_df, alt_mech_df, pol_df, pol_info_df, **kwargs):
        beachball_plot = BeachballPlot(
            mech_df=mech_df, 
            alt_mech_df=alt_mech_df, 
            pol_df=pol_df, 
            pol_info_df=pol_info_df,
            **kwargs)
        
        # Connect the clicked_station_signal to store_event_station
        beachball_plot.clicked_station_signal.connect(
            lambda _ , station_name: self.store_event_station(self.event_id, station_name, self.pick_pol_df, **kwargs))
        
        return beachball_plot
    
    def store_event_station(self, event_id, station_name, pick_pol_df, **kwargs):
        self.event_id = event_id
        self.station_name = station_name
        self.q_update_waveform(self.event_id, self.station_name, self.pick_pol_df, self.mseed_dir, **kwargs)
    
    def q_update_waveform(self, event_id, station_name, pick_pol_df, data_dir, **kwargs):
        self.wf_plot.update_waveform(event_id, station_name, pick_pol_df, data_dir, **kwargs)

    def idx_in_pick_pol_df(self, event_id, station_name, df):
        '''
        *** Note: Pass the original (whole) picks DataFrame, not the one only for the current event
        '''
        return df.loc[(df['event_id'] == event_id) &
                      (df['phase_type'] == 'P') &
                      (df['station_id'].str.split('.').str[1] == station_name.split('.')[0])].index.values[0]

    def reverse_polarity(self):
        if self.event_id and self.station_name:
            # Reverse the polarity in the pick_pol_df DataFrame
            idx = self.idx_in_pick_pol_df(self.event_id, self.station_name, self.pick_pol_df_all)
            print('before', self.pick_pol_df_all.loc[idx, 'phase_polarity'])
            self.pick_pol_df_all.loc[idx, 'phase_polarity'] = -self.pick_pol_df_all.loc[idx, 'phase_polarity']
            self.pick_pol_df.loc[idx, 'phase_polarity'] = -self.pick_pol_df.loc[idx, 'phase_polarity']
            print('after', self.pick_pol_df_all.loc[idx, 'phase_polarity'])

    def set_polarity_zero(self): #, pick_pol_df, event_id, station_name):
        if self.event_id and self.station_name:
            # Set the polarity to zero in the pick_pol_df DataFrame
            idx = self.idx_in_pick_pol_df(self.event_id, self.station_name, self.pick_pol_df_all)
            print('before', self.pick_pol_df_all.loc[idx, 'phase_polarity'])
            self.pick_pol_df_all.loc[idx, 'phase_polarity'] = 0
            self.pick_pol_df.loc[idx, 'phase_polarity'] = 0
            print('after', self.pick_pol_df_all.loc[idx, 'phase_polarity'])

    def save_polarity_file(self):
        reply = QMessageBox.question(
            self, 'Confirmation', 'Are you sure to save the edited polarity file \
                (This will overwrite the existing file)?',
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            self.pick_pol_df_all.to_csv(self.pick_pol_path, index=False)
            print('Polarity file saved successfully')
        
    def rerun_skhash(self, **kwargs):
        reply = QMessageBox.question(
            self, 'Confirmation', 'Are you sure to re-run SkHASH \
                (This will overwrite the existing file)?',
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No)
        if reply == QMessageBox.Yes:
            rerun_skhash = RerunSKHASH(**kwargs)
            # Wait for the process to finish and reload the data
            self.load_data(**kwargs)


# Example usage
if __name__ == '__main__':
    import os 

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
