#To test the code example datasets have been placed in the current directory. 

#The code is used to get the UNIX time of the first gate boundary considering the start time of the correlator is defined and a previous reference profile pbservation has been made. 

#The test is to check the working of the predictor code and verify the predicted phase of the TOA by comparing it with the observed phase at some time in the future. 

#Observation Details:
1.Pulsar J0814+7429 was observed using GMRT on 2 seperate days 31st Jan and 1st Feb 2024. The profiles were genrated using gptool. The hdr files for the respctive observation gives the time for the first phase bin in each profile. 

2. The Pulsar parameters are present in the par file. 

# Test 
1.To predict the phase of the TOA(highest amplitude bin in the profile) given a start time for the correlator.Here we assume that the time in HDR file for observation on 1st Feb is the start time which will help us to compare the predicted phase to the observed phase. 
2.To get the maximum delayed phase of the TOA at highest observing frequency.(dispersion delay)
3. To get the first gate boundary and gate width depending on the pulse width in UNIX time format.  

# Check the parameters in the config file 

If everything is working properly the code should create the following output in the terminal : 

Predicting TOA for Pulsar : J0814+7429 after start time of correlator : 60341.58438817662 

--------Observation Parameters------------
Profile observation time :MJD_observation_UTC:60340.565675949074
Start time of correlator :60341.58438817662
----------Observed TOA from Profile-----------
Observed TOA for J0814+7429 MJD1: 60340.56567637005
Corresponding Fractional phase: 0.0281457994133234
-------Running Tempo2 Polyco-----------
This program comes with ABSOLUTELY NO WARRANTY.
This is free software, and you are welcome to redistribute it
under conditions of GPL license.

[tempo2Util.C:396] Warning: [MISC1] Unknown parameter in par file:  RM
[tempo2Util.C:401] Warning: [DUP1] duplicated warnings have been suppressed.
************************************************* 
Warning: you are running in tempo1 emulation mode 
************************************************* 
WARNING: tzrmjd not set.  Setting to 49162
WARNING: tzrfrq not set.  Setting to 1400
WARNING: tzrsite not set.  Setting to GMRT
[tempo2Util.C:396] Warning: [CLK3] no clock corrections available for clock UTC(GM) for MJD 49162.0
[tempo2Util.C:396] Warning: [CLK4] Trying assuming UTC = UTC(GM)
[tempo2Util.C:396] Warning: [CLK9] ... ok, using stated approximation 
[tempo2Util.C:396] Warning: [CLK6] Proceeding assuming UTC =  UTC(GM)
[tempo2Util.C:396] Warning: [CLK10] MJD is outside of UT1 table range  (MJD 60340.315941)
Created Polyco File :J0814+7429polyco_new.dat
-----------Loading Polyco Coeffs and Rphases------------
DM value: 5.75066
----------Observed TOA from Profile-----------
Observed TOA for J0814+7429 MJD1: 60340.56567637005
Corresponding Fractional phase: 0.0281457994133234
------Final Output-----------
Predicted TOA from start time 2024-02-01 14:01:31.138: 60341.58442145627 at 550.049115913556 MHz
Corresponding Fractional phase: 0.22503459453582764
The predicted pulsar period at 2024-02-01 14:01:31.138 is 1.2922772259871718 sec
---------Gating Output---------
The Toa at highest frequency 750 MHz: 60341.58442103455
The pulse arrival at 750MHz in IST : 2024-02-01 19:31:33.977385
The pulse arrival at 750MHz in UNIX timestamp in seconds : 1706815893.977385
-----------Gates---------------
gate at 10%  of period from TOA at 750 MHz is lower_gate boundary in IST:2024-02-01 19:31:33.758661 in UTC: 60341.58441850302 
gate at 10%  of period from TOA at 750 MHz is higher_gate boundary in IST:2024-02-01 19:31:34.196110 in UTC: 60341.58441850302 
UNIX Timestamp of lower gate boundary in seconds : 1706815893.7586608
UNIX Timestamp of higher gate boundary in seconds : 1706815894.1961095
gate_width from TOA:218.72414468621207 in ms 


# Validity Check:
Plot the profile for 1st Feb using gnuplot and check the phase for the highest amplitude point. 
Verify if the phase is comparable to the value : Corresponding Fractional phase: 0.22503459453582764 

