function conf = input_MARS(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.MARSFile = fullfile(conf.OGCM.BaseDir,'champs_ex.nc');

%%%
