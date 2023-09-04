function conf = input_conf(varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inputs data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf.Data.BaseDir = fullfile(pwd,'InputData'); %Folder containing inputs
conf.Data.Domain = fullfile(conf.Data.BaseDir,'coastline_Gironde.txt'); % domain file (full path)
conf.Data.ParticlesFile = fullfile(conf.Data.BaseDir,'initial_particle_location_4.csv'); % particle release file (full path)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% OGCM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf.OGCM.Model_name = 'TELEMAC'; % Select POM/FVCOM/MARS/TELEMAC/MOHID

conf.OGCM.DomainName = 'Gironde';

conf.OGCM.BaseDir = fullfile(pwd,'InputData');
conf.OGCM.TimeStep =20/60/24; %in days
conf.OGCM.Coordinates = 'cartesian'; %Options: 'spherical' or 'cartesian'
conf.OGCM.VerticalLayer = 'sigma2depthVar'; % Only for 3D mode
                                          % Options: 'rectangular'
                                          %          'sigma2depthVar': sigma layer varing with elevation
                                          %          'hybrid': sigma-rectangular  
conf.OGCM.BottomType = 'Cte'; % Options: 'Cte' for constant bed over time, i.e, no morphological changes
                              %          'Var' for variable bed over time, i.e, morphological changes

% Select an specific region of the grid to increase the resolution (optional, no: default)

conf.OGCM.cut = 'yes'; %'yes' or 'no'

if strcmpi(conf.OGCM.cut,'yes')
%     conf.OGCM.minLon = 3.47*10^5; % min lon
%     conf.OGCM.maxLon = 4.8*10^5; % max lon
%     conf.OGCM.minLat = 6.38*10^6; % min lat
%     conf.OGCM.maxLat = 6.528*10^6; % max lat
    conf.OGCM.minLon = 3.47*10^5; % min lon
    conf.OGCM.maxLon = 4.25*10^5; % max lon
    conf.OGCM.minLat = 6.44*10^6; % min lat
    conf.OGCM.maxLat = 6.528*10^6; % max lat
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Trajectory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mode
conf.Traj.Verbose = 2; % 0: no information, 1: basic informations, 2: full informations
conf.Traj.Mode = '3D'; %'2D' or '3D' tracking
conf.Traj.NbCores = 4 ; % Number of cores for particule parallel processing

% Time parameters
conf.Traj.ReleaseTime = datenum(2021,06,02,0,0,0);
conf.Traj.TrajectoryDuration = 1; % [days]
conf.Traj.TimeStepCalc = 10/60/24; % [days]
%conf.Traj.TimeStepCalc = 5/60/24;
conf.Traj.Direction = 'forward'; %forward or backward
conf.Traj.Scheme = 'RK4'; % Euler: Euler 1st order, RK2 and RK4: Runge Kutta 2nd and 4th order

% Output
conf.Traj.TimeStepOut=30/60/24; % time step outputs
%conf.Traj.TimeStepOut=10/60/24;
conf.Traj.ScenarioName = [conf.OGCM.DomainName '_chunk0.5_TSOut30_TSCalc10_AdvOnly']; % Free choice
conf.Traj.BaseDir = fullfile(pwd,'Outputs'); % Folder to save outputs

% Chunk Method: to avoid memory, the time domain can be divided in
% different partitions (chunks) with a duration of duration of "chunklen" 
% (e.g. 5 days) to calculate trajectories. Necessary for long simulations
% and/or high resolution grids for the 3D mode.
conf.Traj.chunklen = 0.5; % duration of partitions (in days)                         

% Dispersion
conf.Traj.KvOption='Cte'; %Options: 'Cte'=constant; 'fromOGCM'=from hydrodynamic model
conf.Traj.KhOption='Cte'; %Options: 'Cte'=constant; 'fromOGCM'=from hydrodynamic model
if strcmpi(conf.Traj.KhOption,'Cte');
    conf.Traj.Kh = 0.2; % [m2/s]  %Kh=0 no dispersion
end
if strcmpi(conf.Traj.KvOption,'Cte');
    conf.Traj.Kv = 0.00001; % [m2/s] For 3D mode %Kv=0 no dispersion
end

% Other Transport Processes

conf.Traj.Beaching = 'no'; % Beaching: 'yes' or 'no'
conf.Traj.Refloating = 'no'; % Refloating: 'yes' or 'no'
if strcmpi(conf.Traj.Refloating,'yes')
    conf.Traj.Tw = 1; % half-live of particle to remain on the beach before washing off again (in days)
    conf.Traj.RefloatingAtHighTide = 'yes';
    if strcmpi(conf.Traj.RefloatingAtHighTide,'yes')
      conf.Traj.xcoord_elev = -8.93; % x-coord/lon to calculate the time series of elevation
      conf.Traj.ycoord_elev = 42.53; % y-coord/lat to calculate the time seriss of elevation
    end
end
conf.Traj.Deposition = 'no'; % Deposition: 'yes' or 'no'
conf.Traj.Resuspension = 'no'; % Resuspension: 'yes' or 'no'
if strcmpi(conf.Traj.Resuspension,'yes')
  conf.Traj.Sliding = 'no'; % Sliding (bedload): 'yes' or 'no'
  conf.Traj.ResOption = 'waldschlager'; % Options: soulsby, waldschlager,value
    if strcmpi(conf.Traj.ResOption,'value')
        conf.Traj.tau1=NaN;
        conf.Traj.tau2=0;
    end
end


% conf.Traj.Windage = 'no'; % Windage: 'yes' or 'no',  Not available for the moment

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Behaviour (3D mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Water and bed sediment properties (useful for resuspension, sliding and ws calculation)
conf.Beh.WaterDensity = 1.025; % [g cm^-3]  (for resuspension, sliding, and ws calculation)
conf.Beh.WaterViscosity = 10^(-6); % [m^2/s] (for ws calculation and sliding)
conf.Beh.SediD50 = 0.35*10^(-3); % [m] (for resuspension/walschlager and sliding)
conf.Beh.Cm = 0.8;    % For sliding only. Empirical drag coefficient to take into account pressure drag, added mass effect, 
%conf.Beh.Cd = 1;%0.003 ; % For bedload only. Bed drag coefficient

% Particle properties (useful for resuspension, sliding and ws calculation)
conf.Beh.PolymerAcronym ='PET';  % Optional, just for records
conf.Beh.ParticleDensity = 0.90; % [g/cm^3] (for resuspension, sliding, and ws calculation)
% conf.Beh.ParticleDensity = 1.030; % [g/cm^3] (for resuspension, sliding, and ws calculation)
conf.Beh.ParticleSize = 0.001; % Length (largest side: e.g diameter for sphere, length for fibers, longest side for sheets) [m] (for resuspension, bedload, and ws calculation)
conf.Beh.ParticleDequi = 0.0012; % Equivalent diameter (a*b*c)^(1/3) [m] (for resuspension, sliding, and ws calculation)
conf.Beh.ParticleShape='fibre'; % Options: 'sphere','fragment','fibre','sheet' (for ws calculation)

% Settling velocity setup
% Buoyancy (settling/rising)
conf.Beh.InitialWsOption='value'; %Options to define the initial settling velocity:
                                      % 'value': known value, given by use (macro- or microplastic)
                                      % 'formulation': Select a specific formulation (microplastic)
                                      % 'rate': Increasing or decreasing rate
 
if strcmpi(conf.Beh.InitialWsOption,'value') %terminal velocity
%     conf.Beh.Ws0 = 0.0001; % initial settling velocity (0=neutral buoyanci; + value: settling; - value: rising)
    conf.Beh.Ws0 = 0; % initial settling velocity (0=neutral buoyanci; + value: settling; - value: rising)
elseif strcmpi(conf.Beh.InitialWsOption,'formulation')
    conf.Beh.WsForm='zhiyao'; %Options: 'waldschlager': recommended for spheres, fragments and fibers (Waldschlager and Schuttrumpf, 2019; Waldschlager et al., 2019; Jal?n-Rojas et al., 2022)
                                    %         'dellino': recommended for sheets (Jalon-Rojas et al., 2022)
                                    %         'zhiyao'
                                    %         'khatmullina'
    if strcmpi(conf.Beh.WsForm,'waldschlager')
        conf.Beh.CSF=0.00003/((0.005*0.00003)^(1/2)); %Corey Shape Factor: c/(a*b)^(1/2)
        conf.Beh.P=1; %Power roundness
    elseif strcmpi(conf.Beh.WsForm,'khatmullina')
        conf.Beh.Diameter=0.03*10^(-3); %[m]
    elseif strcmpi(conf.Beh.WsForm,'dellino')
        conf.Beh.Phi=0.0324; %Spherecity [-]
    end
elseif strcmpi(conf.Beh.InitialWsOption,'rate')
    conf.Beh.Ws0 = 0.0001;
    conf.Beh.WsRate = 0.005*10^(-3); %mm/day;
end

% Biofouling
conf.Beh.Biofouling='no'; % Options: 'no': no biofouling
                                  %    'inWs': already considered when 
                                  %            defining the ws value
                                  %    'CteCalculated': (stationary) ws increase estimated by TrackMPD
                                  %    'Rate': no stationary, time increase defined by a rate
                                  %    'Empirical': to be added soon
                                        
if strcmpi(conf.Beh.Biofouling,'CteCalculated') % just for spheres or fibers  
    conf.Beh.BiofoulingThickness = 0.01*10^(-3); %meters
    conf.Beh.FilmDensity = 1.15; %g cm^-3
elseif strcmpi(conf.Beh.Biofouling,'Rate')  
    conf.Beh.BiofoulingRate = 0.5*10^(-3); %density rate: g/cm3/day;
end
   

                                        
%Degradation
conf.Beh.Degradation='no'; %Option: 'no': no degradation
                           %         'Rate': size (and ws) decrease define by a rate

if strcmpi(conf.Beh.Degradation,'no')  % just available for 'waldschlager' formulation
     conf.Beh.DegradationRate = -30; %Percentage of size(dequi) decrease per day (negative value)
end                           
                          
 

end
