function conf = input_FVCOM(varargin)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% OGCM inputs, folder and files
conf.OGCM.BaseDir = fullfile(pwd,'InputData'); 
conf.OGCM.FVCOMFile = fullfile(conf.OGCM.BaseDir,'YRE_0001.nc');
conf.OGCM.FVCOMGrid = fullfile(conf.OGCM.BaseDir,'YRE_Geo.grd');
    
% Number of points in the new rectangular grid (lon dimension)    
conf.OGCM.NumLonGrid = 151; 
   
% Number of points in the new rectangular grid (lat dimension)
conf.OGCM.NumLatGrid = 63; 
    
% Name of variable time
conf.OGCM.NameTime = 'time'; 
    
% Name of variable for vertical velocity (3D mode)
conf.OGCM.NameW = 'ww';  

%%%
