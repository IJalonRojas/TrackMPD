% Run TrackMPD v.1

clear all; close all; clc;
conf_name='input_conf'; %Name of input file


%% Step 1
%  Transform the hydrodynamic input files (OGCM outputs) to the “standard” 
%  TrackMPD input format. 

%  OPTIONAL: This transformation may be no necessary if already made in 
%  previous simulations

transformInputs(conf_name)


%% Step 2
%  Run the tracking model

runTrackMPD(conf_name)
