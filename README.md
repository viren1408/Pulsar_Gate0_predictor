# Pulsar Gating Prediction Pipeline

This Python script predicts the Time of Arrival (TOA) and gate boundaries for a pulsar for gated imaging applications. Currently, the code is tailored for GMRT observations 

## Creator

- Viren Mandaogane

### **Caution** 
The code currently is in a developing phase and requires sanity and validity checks which are currently in progress, any suggestions and comments are welcome. 
I am working on a better more general documentation for the code. 


## Description

This script takes observational parameters such as HDR file path, correlator start date and time, and pulsar information including PAR file, pulsar name, and profile paths. It then uses Tempo2 to predict the phase of the pulsar at the given reference time of the correlator. Additionally, it calculates the dispersion delay and estimates the number of rotations of the pulsar between the time of pulsar profile observation and correlator start time. Finally, it outputs the dispersed TOA, predicted TOA, and maximum delay, along with gate boundaries for the pulsar. The gate gets first gate boundaries which can then be used to apply gates in the correlator gate table. 

## Functionality

- **extract_ist_time_date**: Extracts IST time and date from HDR file and converts it to UTC MJD.
- **isttoutc**: Converts IST to UTC time.
- **generate_polyco_tempo**: Generates polycos for a pulsar at a given date and observatory site using Tempo2.
- **load_polyco**: Loads polyco coefficients and required parameters from the polyco file.
- **calculate_polyco_phase**: Calculates the predicted phase and rotation frequency of the pulsar at a given time.
- **load_profile_data**: Loads observed pulsar profile from a file and extracts time and phase at the peak amplitude point.
- **find_max_amplitude_index**: Finds the index and MJD of the peak amplitude point in the pulsar profile.
- **plot_profiles**: Plots the observed pulsar profile.
- **extract_dm_from_par**: Extracts the Dispersion Measure (DM) value from the par file.
- **calculate_delay_table**: Calculates the delay table and maximum delay for a given frequency range and dispersion measure.
- **estimated_number_of_rotations**: Estimates the number of rotations of the pulsar between the time of pulsar profile observation and correlator start time.
- **mjd_to_ist_astropy**: Converts UTC to IST using Astropy.
- **gate_prediction_pipeline**: Runs the pipeline to predict TOA and gate boundaries for a pulsar.
- **read_config**: Reads parameters from a configuration file.

## Usage

1. Prepare a configuration file `gating_config.ini` with appropriate parameters.
2. Run the script `gate_prediction_pipeline.py`.
3. View the output for dispersed TOA, predicted TOA, maximum delay, and gate boundaries.

## Dependencies
- Tempo2 
- Python 3x
- Astropy
- Pandas
- Matplotlib
- NumPy

## Output

The script generates a log file `{pulsar_name}_pred_toa_log.txt` containing the predicted TOA and gate boundaries. Additionally, it may generate plots of observed pulsar profiles if specified.

## Sanity check. 
The current test for validity is also present in the /test dir 
