"""
Predicts TOA and gate boundaries for a pulsar for gated imaging applications

Args:
    hdr_file_path (str): Path to HDR file.
    corr_observation_datetime (str): Correlator observation datetime. (To check use the hdr time for an observation where prediction needs to be verified)
    pulsar_name (str): Name of the pulsar.
    plot_psr_profiles (bool, optional): Whether to plot pulsar profiles. Defaults to True.

Returns:
    tuple: Dispersed TOA, predicted TOA, maximum delay.
"""

# Importing necessary libraries
import math
import numpy as np
import subprocess
import configparser
import sys 
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time, TimeDelta
import pandas as pd
import sys
from astropy.time import Time
from astropy.coordinates import EarthLocation
import astropy.units as u
from datetime import timedelta
import os

####### Read HDR file and convert the IST time for profile start time and convert it to MJD(UTC) ########
def extract_ist_time_date(hdr_file):
    """
    Extracts IST time and date from HDR file and converts it to UTC MJD.

    Args:
        hdr_file (str): Path to HDR file.

    Returns:
        str: UTC MJD in ISO format.
    """
    ist_time = None
    ist_date = None

    with open(hdr_file, 'r') as file:
        for line in file:
            if line.startswith('IST Time:'):
                ist_time = line.split(':', 1)[1].strip()  # Capture IST time 
            elif line.startswith('Date:'):
                date_parts = line.split(':',1)[1].strip().split(':')  # Capture Date 
                ist_date = f'{date_parts[2]}-{int(date_parts[1]):d}-{int(date_parts[0]):d}'  # Convert to ISO format for astropy 

    if ist_time and ist_date:
        ist_datetime_str = f'{ist_date} {ist_time}'
        ist_time = Time(ist_datetime_str, format='iso', scale='utc').iso

    return ist_time


def isttoutc(ist_datetime): 
    """
    Converts IST to UTC time.

    Args:
        ist_datetime (str): IST datetime.

    Returns:
        str: UTC datetime.
    """
    given_date_time = Time(ist_datetime, format='iso', scale='utc')

    offset = 19800.0  # 5 hours and 30 minutes
    
    # Time difference to subtract
    time_difference = TimeDelta(offset, format='sec')

    # Subtract the time difference
    utc_datetime = given_date_time - time_difference
    return utc_datetime


####### Tempo2 to predict the phase of the pulsar at the given reference time of the correlator########
def generate_polyco_tempo(par_file_name, mjd_start, mjd_end, nspan, ncoeffs, maxha, telescope_code, frequency, psr_name):
    """
    Generates polycos for a pulsar at a given date and observatory site.

    Args:
        par_file_name (str): Path to par file.
        mjd_start (float): MJD start.
        mjd_end (float): MJD end.
        nspan (int): N span.
        ncoeffs (int): N coeffs.
        maxha (int): Max ha.
        telescope_code (str): Telescope code. 
        frequency (float): Frequency.
        psr_name (str): Name of the pulsar.
    """
    subprocess.run(f'tempo2 -f "{par_file_name}" -polyco "{mjd_start} {mjd_end} {nspan} {ncoeffs} {maxha} {telescope_code} {frequency}" -polyco_file "{psr_name}" -tempo1', shell=True)


def load_polyco(filename, ncoeffs):
    """
    Loads polyco coefficients and required parameters from the polyco file.

    Args:
        filename (str): Path to polyco file.
        ncoeffs (int): Number of coefficients.

    Returns:
        tuple: MJD values, reference phase values, rotation frequency value, polyco coefficients.
    """
    nlines_each_set = 2 + ncoeffs // 3  # A condition for convenience in dealing with the polyco format 

    with open(filename, "r") as polyco:
        lines = polyco.readlines()

    number_of_sets = len(lines) // nlines_each_set 

    mjd_values = np.zeros(number_of_sets)
    rphase_values = np.zeros(number_of_sets)
    rotation_freq_values = np.zeros(number_of_sets)
    polyco_coefficients = []

    for i in range(number_of_sets):
        mjd_values[i] = float(lines[i * nlines_each_set].split()[3])  # tmid values 
        rphase_values[i] = float(lines[i * nlines_each_set + 1].split()[0])  # reference phase at each tmid value 
        rotation_freq_values[i] = float(lines[i * nlines_each_set + 1].split()[1])
        
        polyco_lines = [line for line in lines[i * nlines_each_set + 2: i * nlines_each_set + nlines_each_set]]
        polyco_coefficients.append([float(coeff) for line in polyco_lines for coeff in line.split()])
    rotation_freq_value = rotation_freq_values[0]  # Reference rotation frequency of pulsar (F0): Taken from the par file and does not update with time 

    return mjd_values, rphase_values, rotation_freq_value, polyco_coefficients


def calculate_polyco_phase(MJDObs, rphase_values, rotation_freq_value, mjd_values, polyco_coefficients, nspan, ncoeffs):
    """
    Calculates the predicted phase and rotation frequency of the pulsar at a given time (MJDObs).

    Args:
        MJDObs (float): Observation time.
        rphase_values (np.ndarray): Reference phase values.
        rotation_freq_value (float): Rotation frequency.
        mjd_values (np.ndarray): MJD values.
        polyco_coefficients (list): List of polyco coefficients.
        nspan (int): N span.
        ncoeffs (int): N coeffs.

    Returns:
        tuple: Frequency, phase.
    """
    dts = []
    for j in range(len(mjd_values)):
        dt = (MJDObs - mjd_values[j]) * 1440  # timespan in minutes 
        dts.append(dt)
    
    nearest_dt = min(abs(dt) for dt in dts) 
    min_index = dts.index(next(dt for dt in dts if abs(dt) == nearest_dt))
   
    if nearest_dt < nspan:  # The chosen tmid should be within the timespan 
        tmid = mjd_values[min_index] 
        if tmid > MJDObs:
            min_index = min_index - 1  # If the nearest dt gives a tmid greater than MJDObs, then select the previous tmid 
            tmid = mjd_values[min_index]
        dt_in = (MJDObs - tmid) * 1440 
        phase = rphase_values[min_index] + (dt_in * 60 * rotation_freq_value)
        frequency = rotation_freq_value
        for i in range(ncoeffs):
            phase += polyco_coefficients[min_index][i] * math.pow(dt_in, i)
            frequency += 1/60*(i * polyco_coefficients[min_index][i] * math.pow(dt_in, i - 1))
        return frequency, phase
    else:
        return None


##### Load Observed Pulsar Profile in a data frame and get the time and phase of the pulsar at peak amplitude point #####
def load_profile_data(file_path, mjd_start):
    """
    Loads observed pulsar profile from a file and extracts time and phase at the peak amplitude point(TOA of the pulsar)

    Args:
        file_path (str): Path to the profile file.
        mjd_start (float): MJD start.

    Returns:
        pd.DataFrame: Loaded profile data.
    """
    df_profile = pd.read_table(file_path, sep='\t', index_col=None)
    df_profile['time_stamp(ms)'] = df_profile.index * 0.08192
    df_profile['timestamp_MJD'] = mjd_start + (df_profile['time_stamp(ms)'] / (1000 * 86400))
    return df_profile


def find_max_amplitude_index(df_profile):
    """
    Finds the index and MJD of the peak amplitude point in the pulsar profile.

    Args:
        df_profile (pd.DataFrame): Profile data.

    Returns:
        tuple: Index of max amplitude, MJD of max amplitude.
    """
    max_amplitude_index = df_profile['value'].idxmax()
    return max_amplitude_index, df_profile['timestamp_MJD'].iloc[max_amplitude_index]  # Defined as TOA for the pulsar 


def plot_profiles(df_profile, observed_phase, toa_observed ,pulsar_name):
    """
    Plots the observed pulsar profile.

    Args:
        df_profile (pd.DataFrame): Profile data.
        observed_phase (float): Observed phase.
        toa_observed (float): Observed TOA.
        pulsar_name (str): Name of the pulsar.
    """
    plt.figure(figsize=(10, 6))
    #plt.plot(df_profile2['#phase'], df_profile2['value'], label='observed phase psr1 mjd2')
    plt.plot(df_profile['#phase'], df_profile['value'], label='observed phase')
    plt.axvline(observed_phase, linestyle='--', color='green',label=f'TOA MJD = {toa_observed}, phase:{observed_phase}')
    
    #plt.axvline(pred_toa_phase, linestyle=':', color='green',label=f'phase at predicted toa({pred_toa}): {pred_toa_phase}')
    
    plt.legend(loc='upper left', fontsize='small', fancybox=True, framealpha=0.3)
    plt.legend(loc='upper left', fontsize='small', fancybox=True, framealpha=0.3)
    plt.xlabel('Phase')
    plt.ylabel('Amplitude Value')
    plt.title(f'{pulsar_name}: Amplitude Value vs Phase')
    plt.legend()
    plt.grid()
    #plt.savefig(f'{pulsar_name}_.png')
    plt.show()
    
    plt.figure(figsize=(10, 6))
    plt.plot(df_profile['time_stamp(ms)'], df_profile['value'], label='observed phase')
    plt.xlabel('Time(ms)')
    plt.ylabel('Amplitude Value')
    plt.title(f'{pulsar_name}: Amplitude Value vs Timestamp(ms)')
    plt.legend()
    plt.grid()
    plt.show()

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##

#### Calculate the dispersion delay and store the maximum delay for the highest frequency ####

def extract_dm_from_par(par_file_path):
    """
    Extracts the DM value from the par file.

    Args:
        par_file_path (str): Path to the par file.

    Returns:
        float: DM value.
    """
    dm_value = None
    with open(par_file_path, 'r') as file:
        for line in file:
            if line.startswith('DM'):
                dm_value = float(line.split()[1])
                break 
    return dm_value


def calculate_delay_table(frequency_low,frequency_high,nchannels,bandwidth, dispersion_measure):
    """
    Calculates the delay table and maximum delay for a given frequency range and dispersion measure.

    Args:
        frequency_low (float): Lowest frequency.
        frequency_high (float): Highest frequency.
        nchannels (int): Number of channels.
        bandwidth (float): Bandwidth.
        dispersion_measure (float): Dispersion measure.

    Returns:
        tuple: Delay table, maximum delay seconds, reference frequency.
    """
    ref_frequency_obs =(frequency_low+(bandwidth/nchannels)/2.0)
    KDM = 4.148808 * 10**3 * dispersion_measure
    delay_table=[]
    freqinterval = bandwidth / nchannels
    frequencies = np.arange(ref_frequency_obs,frequency_high,freqinterval)
    for f in frequencies:
        delay_table.append(((1.0 / ref_frequency_obs**2)-(1.0 / f**2)) * KDM)
    max_delay_seconds = max(delay_table)
    return delay_table,max_delay_seconds,ref_frequency_obs

## !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ##

def estimated_number_of_rotations(MJD_observation_mjd2,observed_toa_mjd1,ref_frequency):
    """
    Estimates the number of rotations of the pulsar between the time of pulsar profile observation and correlator start time

    Args:
        MJD_observation_mjd2 (float): MJD for correlator start time.
        observed_toa_mjd1 (float): Observed TOA.
        ref_frequency (float): Reference frequency.

    Returns:
        int: Estimated number of rotations.
    """
    estimated_rotations = int((MJD_observation_mjd2 - observed_toa_mjd1) * 24 * 60 * 60 *ref_frequency + 2)
    return  estimated_rotations

def mjd_to_ist_astropy(mjd_utc):
    """
    Converts UTC to IST using Astropy.

    Args:
        mjd_utc (float): MJD in UTC.

    Returns:
        datetime: IST datetime.
    """
    # Define the UTC time using Astropy Time
    time_utc = Time(mjd_utc, format='mjd', scale='utc')

    # Define the Earth location for IST

    # Convert to IST
    time_ist = time_utc + timedelta(hours=5, minutes=30)

    return time_ist

def gate_prediction_pipeline(hdr_file_path,corr_observation_datetime, pulsar_name,plot_psr_profiles=True):
    """
    Runs the pipeline to predict TOA and gate boundaries for a pulsar.

    Args:
        hdr_file_path (str): Path to HDR file.
        corr_observation_datetime (str): Correlator observation datetime.
        pulsar_name (str): Name of the pulsar.
        plot_psr_profiles (bool, optional): Whether to plot pulsar profiles. Defaults to True.

    Returns:
        tuple: Dispersed TOA, predicted TOA, maximum delay.
    """
    file = open(f"{pulsar_name}_pred_toa_log.txt", "w")  

# Declare all the relevant time variables 

    # Time of profile observation  
    time_ist_observation = extract_ist_time_date(hdr_file_path) #declare 
    time_utc_observation = isttoutc(time_ist_observation) 

    # Time for the start time of correlator  
    time_utc_corr_observation = isttoutc(corr_observation_datetime) #declare

    observation_date_mjd1 = time_utc_observation 
    observation_date_mjd2 = time_utc_corr_observation

    observation_time_mjd1 = Time(observation_date_mjd1, scale='utc')
    observation_time_mjd2 = Time(observation_date_mjd2, scale='utc')
 
    MJD_observation_mjd1 = observation_time_mjd1.mjd # MJD for profile observation 
    MJD_observation_mjd2 = observation_time_mjd2.mjd #MJD for correlator start time 

    print(f'Predicting TOA for Pulsar : {pulsar_name} after start time of correlator : {MJD_observation_mjd2} ')
    file.write(f'Predicting TOA for Pulsar : {pulsar_name} after start time of correlator : {MJD_observation_mjd2}\n')
    print('')
    file.write('\n')
    print('--------Observation Parameters------------')
    file.write('--------Observation Parameters------------\n')
    
    print(f'Profile observation time :MJD_observation_UTC:{MJD_observation_mjd1}')
    file.write(f'Profile observation time :MJD_observation_UTC:{MJD_observation_mjd1}\n')

    print(f'Start time of correlator :{MJD_observation_mjd2}')
    file.write(f'Start time of correlator:{MJD_observation_mjd2}\n')    

    df_profile_mjd1 = load_profile_data(profile_path_obs1, MJD_observation_mjd1)

    max_amplitude_index_mjd1, observed_toa_mjd1 = find_max_amplitude_index(df_profile_mjd1)

    observed_phase_at_toa_mjd1 = df_profile_mjd1['#phase'][max_amplitude_index_mjd1]

    # condition to plot the pulsar profile and check the observed TOA 
    if plot_psr_profiles:
        plot_profiles(df_profile_mjd1,observed_phase_at_toa_mjd1,observed_toa_mjd1,pulsar_name)

    print('----------Observed TOA from Profile-----------')
    file.write('----------Observed TOA from Profile-----------\n')
    print(f"Observed TOA for {pulsar_name} MJD1: {observed_toa_mjd1}")
    file.write(f"Observed TOA for {pulsar_name} MJD1: {observed_toa_mjd1}\n")
    print(f"Corresponding Fractional phase: { observed_phase_at_toa_mjd1}")
    file.write(f"Corresponding Fractional phase: { observed_phase_at_toa_mjd1}\n")

    print('-------Running Tempo2 Polyco-----------')
    mjd_start_tempo = MJD_observation_mjd1
    mjd_end_tempo = MJD_observation_mjd2+1.0
    nspan = 60 
    ncoeffs = 12
    max_ha = 12
    site = 'GMRT'
    freq_obs_MHz = (frequency_low+(bandwidth/nchannels)/2.0)

    generate_polyco_tempo(par_file,mjd_start_tempo, mjd_end_tempo, nspan, ncoeffs, max_ha, site, freq_obs_MHz, pulsar_name)

    polyco_file = f'{pulsar_name}polyco_new.dat'

    if os.path.exists(polyco_file):
        print(f'Created Polyco File :{polyco_file}')
    else:
        print(f'Could not create Polyco File :{polyco_file} please check the input config file and par file')

    

    print('-----------Loading Polyco Coeffs and Rphases------------')

    mjd_values,rphase_values,rotation_freq_value, polyco_coefficients = load_polyco(polyco_file,ncoeffs)
    
    
    dm_value = extract_dm_from_par(par_file)
    if dm_value is not None:
        print("DM value:", dm_value)
    else:
        print("DM value not found in the par file.")

    delay_table,max_delay_seconds,ref_frequency_obs = calculate_delay_table(frequency_low,frequency_high,nchannels,bandwidth,dm_value)


    print('----------Observed TOA from Profile-----------')
    file.write('----------Observed TOA from Profile-----------\n')
    print(f"Observed TOA for {pulsar_name} MJD1: {observed_toa_mjd1}")
    file.write(f"Observed TOA for {pulsar_name} MJD1: {observed_toa_mjd1}\n")
    print(f"Corresponding Fractional phase: { observed_phase_at_toa_mjd1}")
    file.write(f"Corresponding Fractional phase: { observed_phase_at_toa_mjd1}\n")
    #print(f"Observed TOA for {pulsar_name} MJD2: {observed_toa_mjd2}")
    #print(f"Corresponding Fractional phase: { observed_phase_at_toa_mjd2}")

    rot_freq_mjd2,pred_phase_at_mjd2 = calculate_polyco_phase(MJD_observation_mjd2,rphase_values,rotation_freq_value, mjd_values, polyco_coefficients, nspan, ncoeffs)
    
    estimated_rotations = estimated_number_of_rotations(MJD_observation_mjd2,observed_toa_mjd1,rotation_freq_value)
    #print(estimated_rotations)
    
    absolute_phase_at_toa_mjd = []
    frac_phase_at_toa_mjd = []
    time_elapsed_values = []
    toa_mjd_values = []

    mjd_in = observed_toa_mjd1
    
    phase_init = (pred_phase_at_mjd2)%1.0 #predicted phase at start time 2 
    
    for i in range(estimated_rotations):
        rot_freq, pred_abs_phase =  calculate_polyco_phase(mjd_in,rphase_values,rotation_freq_value, mjd_values, polyco_coefficients, nspan, ncoeffs)
        #freq = polyco_pint.eval_spin_freq(mjd_in)
        absolute_phase_at_toa_mjd.append(pred_abs_phase)
        period = 1 / rot_freq
        phase = pred_abs_phase - phase_init
        frac_phase_toa = (phase % 1.0)
        time_elapsed = (mjd_in - observed_toa_mjd1) * 1440
        toa_mjd_values.append(mjd_in)
        frac_phase_at_toa_mjd.append(frac_phase_toa)
        time_elapsed_values.append(time_elapsed)
        next_toa_mjd = mjd_in + (period / (24 * 60 * 60))
        mjd_in = next_toa_mjd

    pred_toa = toa_mjd_values[-1]
    #pred_toa = pred_toa[0]
    pred_toa_phase = frac_phase_at_toa_mjd[-1]
    #delay_table,max_delay_seconds,ref_frequency_obs = calculate_delay_table(frequency_low,frequency_high,nchannels,bandwidth, dispersion_measure_value)

    dispersed_toa = pred_toa - max_delay_seconds / (24 * 60 * 60)
    # Calculate Unix timestamp from MJD
    mjd_offset = 40587  # MJD offset for 1970-01-01
    unix_timestamp_dispersed_toa = (dispersed_toa - mjd_offset) * 86400  # 86400 seconds in a day
    
    mjd_utc = dispersed_toa  # Replace with your MJD value
    time_ist = mjd_to_ist_astropy(mjd_utc)
    time_object = Time(time_ist, format='mjd')
    datetime_object = time_object.datetime
    #print(datetime_object)
    
    period = 1/rotation_freq_value
    pulse_width_ms_10 = 43 #ms
    gate_width_ms = pulse_width_ms_10+ period*0.1*1000 # width at 10 % of peak (from ATNF) plus 10% of period     #(100*0.08192) # width at 10 % of peak (from ATNF) plus 100 bins 
    gate_0_1 = dispersed_toa - gate_width_ms/(24 * 60 * 60*1000) #converting gatewidth(ms) to MJD and then subtracting from disperesed toa
    unix_timestamp_gate_0_1 = (gate_0_1 - mjd_offset) * 86400

    gate_0_2 = dispersed_toa + gate_width_ms/(24 * 60 * 60*1000)
    unix_timestamp_gate_0_2 = (gate_0_2 - mjd_offset) * 86400
    
    mjd_utc_gate_01 = gate_0_1
    time_ist_gate_01 = mjd_to_ist_astropy(mjd_utc_gate_01)
    time_object_gate_01 = Time(time_ist_gate_01, format='mjd')
    IST_gate01 = time_object_gate_01.datetime
   
    mjd_utc_gate_02 = gate_0_2
    time_ist_gate_02 = mjd_to_ist_astropy(mjd_utc_gate_02)
    time_object_gate_02 = Time(time_ist_gate_02, format='mjd')
    IST_gate02 = time_object_gate_02.datetime

    #formatted_datetime = example_datetime.strftime("%Y-%m-%d %H:%M:%S.%f")
    print('------Final Output-----------')
    file.write('------Final Output-----------\n')
    print(f"Predicted TOA from start time {observation_date_mjd2}: {pred_toa} at {ref_frequency_obs} MHz")
    file.write(f"Predicted TOA from start time {observation_date_mjd2}: {pred_toa} at {ref_frequency_obs} MHz\n")
    print(f"Corresponding Fractional phase: {pred_toa_phase}")
    file.write(f"Corresponding Fractional phase: {pred_toa_phase}\n")
    #print(f"Offset in ms = {(pred_toa-observed_toa_mjd2)*24*60*60*1000}")

    print('---------Gating Output---------')
    file.write('---------Gating Output---------\n')
    print(f'The Toa at highest frequency {frequency_high} MHz: {dispersed_toa}')
    file.write(f'The Toa at highest frequency {frequency_high} MHz: {dispersed_toa}\n')
    print(f'The pulse arrival at {frequency_high}MHz in IST : {datetime_object}') #Date and time : UTC
    file.write(f'The pulse arrival at {frequency_high}MHz in IST : {datetime_object}\n') #Date and time : UTC
    print(f'The pulse arrival at {frequency_high}MHz in UNIX timestamp in seconds : {unix_timestamp_dispersed_toa+19799.99999976158}') #Date and time : UTC
    file.write(f'The pulse arrival at {frequency_high}MHz in UNIX timestamp in seconds : {unix_timestamp_dispersed_toa+19799.99999976158}\n') #Date and time : UTC
    file.write('\n')
    print('-----------Gates---------------')
    file.write('-----------Gates---------------\n')
    print(f'gate at 10%  of period from TOA at 750 MHz is lower_gate boundary in IST:{IST_gate01} in UTC: {gate_0_1} ')
    file.write(f'gate at 10%  of period from TOA at 750 MHz is lower_gate boundary in IST:{IST_gate01} in UTC: {gate_0_1} \n')
    print(f'gate at 10%  of period from TOA at 750 MHz is higher_gate boundary in IST:{IST_gate02} in UTC: {gate_0_1} ')
    file.write(f'gate at 10%  of period from TOA at 750 MHz is higer_gate boundary IST:{IST_gate02} in UTC: {gate_0_2} \n')    
    print(f'UNIX Timestamp of lower gate boundary in seconds : {unix_timestamp_gate_0_1+19799.99999976158}')
    file.write(f'UNIX Timestamp of lower gate boundary in seconds : {unix_timestamp_gate_0_1+19799.99999976158}\n')
    print(f'UNIX Timestamp of higher gate boundary in seconds : {unix_timestamp_gate_0_2+19799.99999976158}')
    file.write(f'UNIX Timestamp of higher gate boundary in seconds : {unix_timestamp_gate_0_2+19799.99999976158}\n')
    print(f'gate_width from TOA:{gate_width_ms}')
    file.write(f'gate_width from TOA:{gate_width_ms} in ms \n')

    #plot_profiles(df_profile_mjd1, observed_phase_at_toa_mjd1,observed_toa_mjd1, pred_toa, pred_toa_phase, pulsar_name)
    #plot_profiles(df_profile_mjd1, df_profile_mjd2, observed_phase_at_toa_mjd1, observed_phase_at_toa_mjd2,observed_toa_mjd1, observed_toa_mjd2, pred_toa, pred_toa_phase, pulsar_name)

    file.close()
    
    return dispersed_toa , pred_toa ,max_delay_seconds

import configparser

def read_config(file_path):
    config = configparser.ConfigParser()
    config.read(file_path)
    return config

config = read_config('gating_config.ini')

# Extracting observational parameters
psr_hdr_file_path = config['observation_dates']['hdr_file_path']
correlator_start_date_time = config['observation_dates']['observation_date2']
par_file = config['pulsar_info']['par_file']
pulsar_name = config['pulsar_info']['pulsar_name']
profile_path_obs1 = config['pulsar_info']['profile_path_obs1']
frequency_low = int(config['pulsar_info']['frequency_low'])
frequency_high = int(config['pulsar_info']['frequency_high'])
bandwidth = int(config['pulsar_info']['bandwidth'])
nchannels = int(config['pulsar_info']['nchannels'])
# Call the function with the extracted parameters
dispersed_toa_psr1, pred_toa_psr1, max_delay_seconds_psr1 = gate_prediction_pipeline(psr_hdr_file_path, correlator_start_date_time, pulsar_name, plot_psr_profiles=True)

