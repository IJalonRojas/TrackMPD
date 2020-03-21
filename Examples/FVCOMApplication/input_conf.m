function conf = inputs_conf(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf.Data.BaseDir = [pwd '\InputDataFVCOM'];
conf.Data.Domain = [conf.Data.BaseDir '\coastline_YRE.dat'] ;
conf.Data.ParticlesFile = [conf.Data.BaseDir '\Initial_particle_locationsYRE.csv']; %first data:particles


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf.OGCM.Model_name = 'FVCOM'; %Select POM or FVCOM

conf.OGCM.DomainName = 'YRE';
conf.OGCM.BaseDir = [pwd '\InputDataFVCOM'];
conf.OGCM.TimeStep = 1/24; %in days
conf.OGCM.VerticalLayer = 'sigma2depthCte'; % Only for 3D mode
                                          % Options: 'Rectangular'
                                          %          'Sigma2depthCte': cte over time, independient of Elevation
                                          %          'Sigma2depthVar': variable over time,dependient of Elevation

if strcmpi(conf.OGCM.Model_name,'POM') % When using transformPOMoutputs function
    
    %%%Particular inputs for POM (SARCCM version)
    conf.OGCM.POM_Prefix = 'u_h';
    conf.OGCM.POM_Suffix = '.nc';
    conf.OGCM.t0 = datenum(2018,01,01,00,0,0); %POM reference time variable (first time step)
    conf.OGCM.Hmin = 1; %H for land
    %%%
    
elseif strcmpi(conf.OGCM.Model_name,'FVCOM') % When using transformFVCOMoutputs function
    
    %%%Particular inputs for FVCOM (SARCCM version)
    conf.OGCM.FVCOMFile = [conf.OGCM.BaseDir '\YRE_0001.nc'];
    conf.OGCM.FVCOMGrid = [conf.OGCM.BaseDir '\YRE_Geo.grd'];
    conf.OGCM.NumLonGrid = 151; %number of points in the new rectangular grid (lon dimension)
    conf.OGCM.NumLatGrid = 63; %number of points in the new rectangular grid (lat dimension)
    conf.OGCM.NameTime = 'time'; %name of variable time
    conf.OGCM.NameW = 'ww';  %name of variable for vertical velocity (3D mode)
    %%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Mode
conf.Traj.Mode = '3D'; %'2D' or '3D' tracking

% Time parameters
conf.Traj.ReleaseTime = datenum(2011,05,21,2,0,0);
conf.Traj.TrajectoryDuration = 1; % in days
conf.Traj.TimeStep = 20/24/60; % in days
conf.Traj.Direction = 'forward'; %forward or backward

% Output
conf.Traj.ScenarioName = [conf.OGCM.DomainName '_example_estuary']; %Free choice
conf.Traj.BaseDir = [pwd '\Outputs']; %Folder to save outputs

% Chunk Method: to avoid memory, the time domain can be divided in
% different partitions (chunks) with a duration of duration of "chunklen" 
% (e.g. 5 days) to calculate advection. Necessary for long simulations
% and/or high resolution grids for the 3D mode.
conf.Traj.chunklen = 1; % duration of partitions (in days)                      

% Dispersion
conf.Traj.Kh = 1; %m2/s
conf.Traj.Kv = 0.00001; %0.00001; %m2/s For 3D mode

% Beaching
conf.Traj.Beaching = 'yes'; %'yes' or 'no'

% Refloating
conf.Traj.Refloating = 'no'; %'yes' or 'no'
if strcmpi(conf.Traj.Refloating,'yes')
    conf.Traj.Tw = 2; % half-live of particle to remain on the beach before washing off again (in days)
    conf.Traj.xcoord_elev = 80; %grid cell (x-coord) to calculate the time series of elevation
    conf.Traj.ycoord_elev = 28; %grid cell (y-coord) to calculate the time seriss of elevation
end

% Windage
conf.Traj.Windage = 'no'; %'yes' or 'no'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Behaviour (3D mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %BehaviourType:
    % 1: macro or microplastic, plastic density<1 g cm^-3, no biofouling, no degradation
    % 2: macroplastic (or microplastic with a known ws), plastic density>1 g cm^-3, no biofouling, no degradation
    % 3: microplastic (shape,size,density), plastic density>1 g cm^-3, no biofouling, no degradation
    % 4: microplastic (shape,size,density), Biofouling (stationary), no degradation
    % 5: microplastic (shape,size,density), Biofouling (non-stationary), no degradation
    % 6: microplastic (shape,size,density), No Biofouling, degradation

BehaviourType = 1;
 
if BehaviourType==1
    conf.beh=behaviour1;
    
    %%%Inputs (Useful if windage)
    conf.beh.PolymerDensity = 0.9; %g cm^-3
    %%%
    
elseif BehaviourType==2
    conf.beh=behaviour2;
    
    %%%Inputs
    conf.beh.Ws = 0.002;
    %%%
    
    
elseif BehaviourType==3
    conf.beh=behaviour3;
    
    %%%Inputs
    conf.beh.PolymerType = 'acrylonitrile-butadiene-styrene'; %Not mandatory
    conf.beh.PolymerAcronym ='ABS'; %Not mandatory
    conf.beh.PolymerDensity = 1.665; %g cm^-3
    conf.beh.Category = 'LongCylinder'; % (Sphere, LongCylinder,ShortCylinder)
    conf.beh.Shape = Cylinder(0.15*10^(-3),0.3*10^(-3)); %Sphere(R) or Cylinder(R,L) (m)
    
    conf.beh.WaterDensity = 1.025; %g cm^-3)
    conf.beh.WaterViscosity = 10^(-6); % m^2/s
    %%%
    
    
elseif BehaviourType==4
    conf.beh=behaviour4;
    
    %%% Inputs
    conf.beh.PolymerType = 'acrylonitrile-butadiene-styrene'; %Not mandatory
    conf.beh.PolymerAcronym = 'ABS'; %Not mandatory
    conf.beh.PolymerDensity = 1.026; %g cm^-3
    conf.beh.Category = 'Sphere'; % (Sphere, LongCylinder,ShortCylinder)
    conf.beh.Shape = Sphere(0.15*10^(-3)); %Sphere(R) or Cylinder(R,L) (m)
    conf.beh.BiofoulingThickness = 0.01*10^(-3); %meters
    conf.beh.FilmDensity = 1.05; %g cm^-3
    
    conf.beh.WaterDensity = 1.025; %g cm^-3
    conf.beh.WaterViscosity = 10^(-6); % m^2/s
    %%%

    
elseif BehaviourType==5
    
    conf.beh=behaviour5;
    
    %%% Inputs
    conf.beh.PolymerType = 'acrylonitrile-butadiene-styrene'; %Not mandatory
    conf.beh.PolymerAcronym = 'ABS'; %Not mandatory
    conf.beh.PolymerDensity = 1.026; %g cm^-3
    conf.beh.Category = 'LongCylinder'; % (Sphere, LongCylinder,ShortCylinder)
    conf.beh.Shape = Cylinder(0.15*10^(-3),0.3*10^(-3)); %Sphere(R) or Cylinder(R,L) (m)
    
    conf.beh.BiofoulingThickness0 = 0*10^(-3); %meters
    conf.beh.FilmDensity = 1.05; %g cm^-3
    conf.beh.BiofoulingRate = 0.005*10^(-3); %meters/day;
    
    conf.beh.WaterDensity = 1.025; %g cm^-3
    conf.beh.WaterViscosity = 10^(-6); % m^2/s
    %%%
   
   
    
elseif BehaviourType==6
     
    conf.beh=behaviour6;
    
    %%%Inputs
    conf.beh.PolymerType = 'acrylonitrile-butadiene-styrene'; %Not mandatory
    conf.beh.PolymerAcronym = 'ABS'; %Not mandatory
    conf.beh.PolymerDensity = 1.026; %g cm^-3
    conf.beh.Category = 'Sphere'; % (Sphere, LongCylinder,ShortCylinder)
    conf.beh.Shape0 = Sphere(1*10^(-3)); %Sphere(R) or Cylinder(R,L) (m)
    conf.beh.DegradationRate = 30; %Percentage of size decrease per day
    
    conf.beh.WaterDensity = 1.025; %g cm^-3)
    conf.beh.WaterViscosity = 10^(-6); % m^2/s
    %%%
    

end
 

end