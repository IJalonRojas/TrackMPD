% Run TrackMPD v.2

clear all; 
conf_name='input_conf'; % Name of input file
confOGCM_name='input_TELEMAC'; 

%% Step 1
%  Transform the hydrodynamic input files (OGCM outputs) to the standard TrackMPD input format. 
%  OPTIONAL: This transformation may be not necessary if already made in previous simulations
tic
transformInputs(conf_name,confOGCM_name)
toc

%% Step 2
%  Run the tracking model
tic
runTrackMPD(conf_name)
toc

