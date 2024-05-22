function runTrackMPD(conf_name,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script run TrackMPD model (2D and 3D mode)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: runTrackMPD_3D.m 001 2018-02-12 10:55:10Z ef $
%        Last version: 2023-03-10 vmarieu
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas & Vincent Marieu 
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('running TrackMPD...')

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic info on configuration of the hydrodynamic model inputs,
% trajectory, particles, domain and grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Configuration
conf=feval(conf_name);

% Trajectory info
ReleaseTime=conf.Traj.ReleaseTime; % Initial drifters release date

switch lower(conf.Traj.Direction) % direction
  case 'backward'
    direction=-1;
  case 'forward'
    direction=1;
  otherwise
    error('wrong direction in inputs_conf')
end

t0=ReleaseTime; %start date
trajDuration=conf.Traj.TrajectoryDuration;
tend=ReleaseTime+direction*trajDuration; %final date


% Particles
pX  =  1;               % Particle X-coordinate
pY  =  2;               % Particle Y-coordinate
pZ  =  3;               % Particle Z-coordinate

parfile = conf.Data.ParticlesFile;

M = dlmread(parfile,',',0,0);
par(:,pX)=M(:,1); % Create variable "par" with particles information
par(:,pY)=M(:,2);
par(:,pZ)=M(:,3);

numpar=size(par,1);
disp(['Number of particles = ' num2str(numpar) ' ; Release Time = ' datestr(ReleaseTime,'dd-mmm-yyyy HH:MM') ' ; ' conf.Traj.Direction ' trajectory']);

% Parallel pool management
if conf.Traj.NbCores > numpar
  disp(['Warning : conf.Traj.NbCores > Number of particles, conf.Traj.NbCores is reduced accordingly:'])
  conf.Traj.NbCores = numpar;
end
if mod(numpar,conf.Traj.NbCores) > 0
  disp(['Warning : The number of particles is not a multiple of conf.Traj.NbCores, it is not optimal !'])
end
if conf.Traj.NbCores > maxNumCompThreads('automatic')
  disp(['Warning : conf.Traj.NbCores > Number of cores on the machine (' num2str(NbCores) '), NbCores is reduced accordingly:'])
  conf.Traj.NbCores = maxNumCompThreads('automatic');
end
disp(['conf.Traj.NbCores: ' num2str(conf.Traj.NbCores)])
Pool = gcp('nocreate');
if isempty(Pool)
  parpool('local',conf.Traj.NbCores);
else
  if Pool.NumWorkers~=conf.Traj.NbCores
    delete(Pool);
    parpool('local',conf.Traj.NbCores);
  end
end

% Test on time steps
if conf.Traj.TimeStepCalc > conf.Traj.TimeStepOut
  error('conf.Traj.TimeStepCalc > conf.Traj.TimeStepOut !')
elseif mod(conf.Traj.TimeStepOut,conf.Traj.TimeStepCalc)~=0
  disp(['Warning : conf.Traj.TimeStepOut is not a multiple of conf.Traj.TimeStepCalc, conf.Traj.TimeStepCalc is adjusted:'])
  TimeStepFactor = ceil(conf.Traj.TimeStepOut/conf.Traj.TimeStepCalc);
  disp(['conf.Traj.TimeStepCalc= ' num2str(conf.Traj.TimeStepOut/TimeStepFactor) ' instead of ' num2str(conf.Traj.TimeStepCalc)])
  conf.Traj.TimeStepCalc = conf.Traj.TimeStepOut/TimeStepFactor;
end

% Particule behaviour configuration
if strcmpi(conf.Traj.Refloating,'yes') 
  if strcmpi(conf.Traj.Beaching,'no')
    error('Refloating cannot occurr without beaching !')
  end
end
if strcmpi(conf.Traj.Resuspension,'yes')
  if strcmpi(conf.Traj.Deposition,'no')
    error('Resuspension cannot occurr without deposition !')
  elseif strcmpi(conf.Traj.ResOption,'soulsby')
    taucr1= Soulsby_tcr1(conf.Beh.ParticleSize,conf.Beh.ParticleDensity,conf.Beh.WaterDensity); %N/m2
    taucr2= Soulsby_tcr2(conf.Beh.ParticleSize,conf.Beh.ParticleDensity,conf.Beh.WaterDensity); %N/m2
  elseif strcmpi(conf.Traj.ResOption,'waldschlager')
    taucr1_= Soulsby_tcr1(conf.Beh.SediD50,2.65,conf.Beh.WaterDensity); %N/m2
    taucr2_= Soulsby_tcr2(conf.Beh.SediD50,2.65,conf.Beh.WaterDensity); %N/m2
    taucr1 = exp_hid(taucr1_,conf.Beh.ParticleDequi,conf.Beh.SediD50);
    taucr2 = exp_hid(taucr2_,conf.Beh.ParticleDequi,conf.Beh.SediD50);
  elseif strcmpi(conf.Traj.ResOption,'value')
      taucr1 = conf.Traj.tauc1;
      taucr2 = conf.Traj.tauc2;
  end
else
  taucr1 = 0.0; % Necessary for parallel computation even if not used
  taucr2 = 0.0;
end

% Domain
if conf.Traj.Verbose >= 2
  disp('Domain and grid loading...')
end
domain=load(conf.Data.Domain);
xland=domain(:,1);
yland=domain(:,2);

% Hydrodynamic model grid
infogrid=load(fullfile(conf.Data.BaseDir,'grid.mat'));
longitude=infogrid.Lon; 
latitude=infogrid.Lat; 
H=infogrid.BottomDepth;
[Lon,Lat]=meshgrid(longitude,latitude); 
mask_water=infogrid.mask_water;
DLon = min(longitude(2:end)-longitude(1:end-1));
DLat = min(latitude(2:end)-latitude(1:end-1));
CoordType=conf.OGCM.Coordinates;
LonDomain = [min(longitude(:)) min(longitude(:)) max(longitude(:)) max(longitude(:))];
LatDomain = [min(latitude(:)) max(latitude(:)) max(latitude(:)) min(latitude(:))];


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define time span of the hydrodynamic and trajectory models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time span of the hydrodynamic inputs
if conf.Traj.Verbose >= 2
  disp('Defining time span...')
end
TimeStepInp=conf.OGCM.TimeStep;
load(fullfile(conf.Data.BaseDir,'timestamps.mat')); % load timestamps
tspanhydro=timestamps;

% time span of the trajectories
TimeStepCalc=conf.Traj.TimeStepCalc;
TimeStepOut=conf.Traj.TimeStepOut;
tspanTOT=t0:direction*TimeStepOut:tend;

% Check if TrackMPD time step is lower than the model time step
TSRatio = ceil(TimeStepOut/TimeStepInp);
if TSRatio > 1
  disp(['Warning: Time Step ratio = ' num2str(TSRatio) ])
  disp('TrackMPD time step (conf.Traj.TimeStepOut) is larger than the model time step (conf.OGCM.TimeStep) !')
  disp('Consider changing to a smaller TrackMPD time step for better performances')
end

% check if the trajectory period (spanTOT) is covered by tspanhydro
% IMPORTANT(***): The calculation of advection in the internal loop needs
% data from two aditional time steps (TimeStepOut) after tend

if direction==1 %forward
  
  if t0<tspanhydro(1) || tend+TimeStepOut>tspanhydro(end) %(***) 1 more time steps
    error ('the hydrodynamic simulation does not cover the particles trajectory')
    %     else
    %         %create a variable with the name of inputs files
    %         for i=1:length(tspanhydro) %forward
    %             OGCMFileNames(i)={[conf.OGCM.BaseDir '\TrackMPDInput' num2str(i) '.mat']};
    %         end
    
  end
  
elseif direction==-1 %backward
  
  if t0>tspanhydro(end) || tend-TimeStepOut<tspanhydro(1) %(***) 1 more time steps
    error ('the hydrodynamic simulation does not cover the particles trajectory')
    %     else
    %         %create a variable with the name of inputs files
    %         for i=1:length(tspanhydro)%backward
    %             OGCMFileNames(i)={[conf.OGCM.BaseDir '\TrackMPDInput' num2str(length(tspanhydro)-i+1) '.mat']};
    %         end
  end
end




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Chunk Method Parameters(IJR 23/05/2018)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%chunklen=conf.Traj.chunklen; %Number of inputs files to be read per chunk
%chunklen_time=chunklen*TimeStepInp; %chunk lengh duration in days

chunklen_time=conf.Traj.chunklen; %chunk lengh duration in days

%Checks: the chunk lenght can't be longer than the trajectory period
%        chunks with only one value are also avoided
%assert(~rem(trajDuration-TimeStepOut,chunklen_time),'bad')

if chunklen_time>trajDuration
  chunklen_time=trajDuration;
  disp('Warning : chunklen exceeds the trajectory duration')
  disp('Warning : chunklen was adjusted to the trajectory duration')
elseif rem(trajDuration-TimeStepOut,chunklen_time)==0
  chunklen_time=chunklen_time+TimeStepOut;
  disp('Warning : chunklen was modified to chunlen+1 in order to meet rem(NTOTFiles-1,chunklen)~=0')
end


ntimesteps_chunk=round(chunklen_time/TimeStepOut)+1; %number of TimeStepOut per chunk

ichunk=1:ntimesteps_chunk:length(tspanTOT);
if mod(length(tspanTOT),ntimesteps_chunk)==0
  jchunk=ntimesteps_chunk:ntimesteps_chunk:length(tspanTOT); %final chunks intervals
else
  jchunk=[ntimesteps_chunk:ntimesteps_chunk:length(tspanTOT) length(tspanTOT)]; %final chunks intervals
end


IniTimeChunk=tspanTOT(ichunk);
FinalTimeChunk=tspanTOT(jchunk);

numpartition = length(ichunk);


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TRAJECTORY CALCULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

% Initialize particles variables
FateType=zeros(numpar,1); %Initial region 0=water; 1=land; 2=bottom; 3=out of Domain
LastUpInput=zeros(numpar,1);
LastVpInput=zeros(numpar,1);
TimeLand=NaN(numpar,1);
TimeOutDomain=NaN(numpar,1); %IJR 14/02/24
TimeSettling=NaN(numpar,1); %IJR 14/02/24
PosWater1=NaN(numpar,1);
PosWater2=NaN(numpar,1);
PosWater3=NaN(numpar,1);
PosBottom1=NaN(numpar,1);
PosBottom2=NaN(numpar,1);
PosBottom3=NaN(numpar,1);

% % Windage: read parameters and velocities
%
% if strcmpi(conf.Traj.Windage,'yes')
%     AirDensity=1.2*10^(-3);
%     WaterDensity=1.025;
%    [DensityRate,SurfaceRate]=windageParameters(AirDensity,WaterDensity,Behaviour.PolymerDensity);
%
%     % IMPORTANT NOTE
%     % We need Uwind Vwind matrix [lat*lon*t] (m/s)
% end

% Start external loop: chunk method
fprintf('Starting chunk method; number of partitions: %i\n',numpartition);
chunk=0; %Initialize chunk counter
%cont=1; %Initialize counter for each time step

for k=1:numpartition % EXTERNAL LOOP
  
  chunk=chunk+1;  %Chunk numeration follows the forward or backward direction
  
  % Select Hydro Input Files for the current chunk
  
  t0chunk=IniTimeChunk(k); %Initial time of chunk
  tendchunk=FinalTimeChunk(k); %Final time of chunk
  tspan_chunk=t0chunk:TimeStepOut*direction:tendchunk; %tspan of chunk
  
  if direction==1 %forward
    
    if chunk==1
      posFile1=find(tspanhydro<=t0chunk,1,'last'); %For the first chunk, this value is not needed
    else
      posFile1=find(tspanhydro<=t0chunk-TimeStepOut,1,'last'); % (***) the ode need 2 more time steps
    end
    posFileEnd=find(tspanhydro>=tendchunk+TimeStepOut,1,'first'); % (***) the ode need 2 more time steps
    
  elseif direction==-1 %backward
    
    if chunk==1
      posFile1=find(tspanhydro>=t0chunk,1,'first');
    else
      posFile1=find(tspanhydro>=t0chunk+TimeStepOut,1,'first');
    end
    posFileEnd=find(tspanhydro<=tendchunk-TimeStepOut,1,'last'); % (***) the ode need 2 more time steps
    
  end
  
  
  NumFilesChunk=abs(posFileEnd-posFile1)+1;
  if conf.Traj.Verbose >= 2
    disp(['Number of files in the chunk: ' num2str(NumFilesChunk)])
  end
  OGCMFileNames=cell(NumFilesChunk,1);
  contFile=1;
  for i=posFile1:direction:posFileEnd
    OGCMFileNames(contFile)={fullfile(conf.Data.BaseDir,['TrackMPDInput' num2str(i) '.mat'])};
    contFile=contFile+1;
  end

  
  % Read OGCM history for the chunk k
  if conf.Traj.Verbose >= 2
    disp('File concatenation...')
  end
  TT=[];
  U=[];
  V=[];
  W=[];
  Depth=[];
  E=[];
  KV=[]; %IJR April 2020
  KH=[]; %IJR April 2020
  for i=1:NumFilesChunk
    
    dataOGCM_aux=load(OGCMFileNames{i});
    TT=cat(1,TT,dataOGCM_aux.time);
    U=cat(4,U,dataOGCM_aux.u);
    V=cat(4,V,dataOGCM_aux.v);
    W=cat(4,W,dataOGCM_aux.w);
    if strcmpi(conf.Traj.KhOption,'fromOGCM')
      KH=cat(4,KH,dataOGCM_aux.kh); %IJR April 2020
    end
    if strcmpi(conf.Traj.KvOption,'fromOGCM')
      KV=cat(4,KV,dataOGCM_aux.kv); %IJR April 2020
    end
    
    if strcmpi(conf.OGCM.VerticalLayer,'sigma2depthVar') || strcmpi(conf.OGCM.VerticalLayer,'hybrid')
      Depth=cat(4,Depth,dataOGCM_aux.depth);
    elseif strcmpi(conf.OGCM.VerticalLayer,'rectangular')
      Depth=cat(4,Depth,dataOGCM_aux.depth);
    else
      error ('Invalid type of vertical layer')
    end
    
    if strcmpi(conf.Traj.Refloating,'yes')==1 || strcmpi(conf.OGCM.VerticalLayer,'hybrid')==1
      % We will need E for hightide option in refloating and to calculate Depth_partLonLat for hybrid vertical layers
      E=cat(3,E,dataOGCM_aux.E);
    end
    
  end
  
  % Velocity unit conversion from m/s to grid units
  if conf.Traj.Verbose >= 2
    disp('Velocity units conversion...')
  end
  [U, V, W] = transformU(Lon, Lat, U, V, W, conf.OGCM.Coordinates); %IJR22
  TheDate = now;
  
  fprintf('******* Current time: %s  === Processing release time: %s Chunk: %s \n', ...
    datestr(TheDate,0), datestr(ReleaseTime,0), num2str(chunk));
  
  % Initialize Refloating variables
  if conf.Traj.Verbose >= 2
    disp('Refloating variables initialization...')
  end
  peak_pos = [];
  if strcmpi(conf.Traj.Refloating,'yes')
    if strcmpi(conf.Traj.RefloatingAtHighTide,'yes')
      Elevation_ts=nan(size(E,3),1);
      for lv=1:size(E,3)
        Elevation_ts(lv)=griddata(Lon,Lat,E(:,:,lv),conf.Traj.xcoord_elev,conf.Traj.ycoord_elev);
      end
      if isnan(sum(Elevation_ts))
        error('conf.Traj.xcoord_elev conf.Traj.ycoord_elev or is not in the domain !')
      end
      Elevation_TRAJ=interp1(TT,Elevation_ts,tspan_chunk); % Interpolate the elevation time series at external time
      
      % Identify high tide
      %[~,peak_pos]=findpeaks(Elevation_TRAJ); % uses signal processing toolbox
      kTime = 0;
      for iTime=2:length(Elevation_TRAJ)-1
        if  Elevation_TRAJ(iTime) > Elevation_TRAJ(iTime-1) && Elevation_TRAJ(iTime) >= Elevation_TRAJ(iTime+1)
          kTime = kTime+1;
          peak_pos(kTime) = iTime;
        end
      end
    end
  end
  
  % Variables for resuspension
  if strcmpi(conf.Traj.Resuspension,'yes')==1
    if conf.Traj.Verbose >= 2
      disp('Resuspension variables initialization...')
    end  
   % difDepth=abs(squeeze(Depth(:,:,end,:)-Depth(:,:,end-1,:)));
   % Ush = squeeze(U(:,:,end-1,:)); % Shear velocity (Lon/day)
   % Vsh = squeeze(U(:,:,end-1,:)); % Shear velocity (Lat/day)
        
    % Bottom layer searched to find Ush
    [Nx,Ny,Nz,Nt]=size(U);
    NearBottomLayer = ones(Nx,Ny);
    Velocity = squeeze(mean(abs(sqrt(U(:,:,:,:).^2 + V(:,:,:,:).^2)),4));
    for iz=1:Nz-1
      NearBottomLayer(Velocity(:,:,iz)>0.0) = iz;
    end
%    figure;
%    pcolor(Velocity(:,:,1));shading flat;colorbar;drawnow
%    figure;
%    pcolor(Velocity(:,:,Nz-1));shading flat;colorbar;drawnow
%    figure;
%    pcolor(NearBottomLayer);shading flat;colorbar;drawnow
    % Shear velocity 
    Ush = zeros(Nx,Ny,Nt);
    Vsh = zeros(Nx,Ny,Nt);
    difDepth=zeros(Nx,Ny,Nt); %IJR 10/10/22
    for ix = 1:Nx
      for iy = 1:Ny
        iz = NearBottomLayer(ix,iy);
        Ush(ix,iy,:) = U(ix,iy,iz,:); 
        Vsh(ix,iy,:) = V(ix,iy,iz,:); 
        difDepth(ix,iy,:)=abs(squeeze(Depth(ix,iy,iz+1,:)-Depth(ix,iy,iz,:)));
      end
    end
    if strcmpi(conf.Traj.KvOption,'fromOGCM')
      Kvsh=squeeze(KV(:,:,end-1,:)); % Kv at last layer  => Voir s'il faut aussi aller chercher à Layer
    else
      Kvsh=conf.Traj.Kv; % Kv at last layer   
    end
  else % Must be set for parallel computing even if not used
    Ush = 0.0;
    Vsh = 0.0;
    difDepth = 0.0;
    Kvsh = 0.0;
  end
  %figure;
  %pcolor(longitude,latitude,sqrt(Ush(:,:,5).^2+Vsh(:,:,5).^2));shading flat;drawnow
    
  %disp(['After preparation: ' num2str(toc,'%.1f') ' s'])
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% INTERNAL LOOP (IJR 25/05/2018
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Initialize variables
  if conf.Traj.Verbose >= 2
    disp('Trajectory variables initialization...')
  end
  len_tspan_chunk=length(tspan_chunk);
  ln_Tot_Output = nan(numpar,len_tspan_chunk);
  lt_Tot_Output = nan(numpar,len_tspan_chunk);
  h_Tot_Output = nan(numpar,len_tspan_chunk);
  LastUpOutput = zeros(numpar,1);
  LastVpOutput = zeros(numpar,1);
  ts_Tot_Output = nan(numpar,len_tspan_chunk); %IJR22
  DepthBottomTraj = nan(numpar,len_tspan_chunk);
  Duration = zeros(1,numpar);
  
  % First value = Release position
  if chunk==1
    ts_Tot_Output(:,1)=ReleaseTime;
    ln_Tot_Output(:,1)=par(:,1);
    lt_Tot_Output(:,1)=par(:,2);
    h_Tot_Output(:,1)=par(:,3);
    
    % For parallel computing, the input arrays cannot be also the output arrays
    ln_Tot_Input = par(:,1);
    lt_Tot_Input = par(:,2);
    h_Tot_Input = par(:,3);
    DepthBottomTraj(:,1)=griddata(Lon(mask_water==1),Lat(mask_water==1),H(mask_water==1),ln_Tot_Output(:,1),lt_Tot_Output(:,1),'nearest'); %Bottom Depth at the release point
    IniLoop=2;
  else
    IniLoop=1;
  end
 
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Internal loop: calculate position for each time step
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if conf.Traj.Verbose >= 2
    disp('Begining parallel loop on particles...')
  end
  parfor (part=1:numpar,conf.Traj.NbCores)
    tic           
  
    ln_Tot_Part = zeros(1,len_tspan_chunk);
    lt_Tot_Part = zeros(1,len_tspan_chunk);
    h_Tot_Part = zeros(1,len_tspan_chunk);
    ts_Tot_Part = zeros(1,len_tspan_chunk);
    DepthBottomTrajPart = zeros(1,len_tspan_chunk);
    
    LL1 = ln_Tot_Input(part);
    LL2 = lt_Tot_Input(part);
    LL3 = h_Tot_Input(part);
    if IniLoop == 2
      ln_Tot_Part(1) = LL1;
      lt_Tot_Part(1) = LL2;
      h_Tot_Part(1) = LL3;
      ts_Tot_Part(1) = ReleaseTime;
    end
    FateTypePart = FateType(part);
    LastUpPart = LastUpInput(part);
    LastVpPart = LastVpInput(part);
    TimeLandPart = TimeLand(part);
    TimeOutDomainPart = TimeOutDomain(part);
    TimeSettlingPart = TimeSettling(part);
    PosWater1Part = PosWater1(part);
    PosWater2Part = PosWater2(part);
    PosWater3Part = PosWater3(part);
    PosBottom1Part = PosBottom1(part);
    PosBottom2Part = PosBottom2(part);
    PosBottom3Part = PosBottom3(part);
        
    % Initialization only to avoid warnings during parallel execution
    ln_ADV = 0;
    lt_ADV = 0;
    h_ADV = 0;
    Up = 0;
    Vp = 0;
    Depth_partLonLat = 0.0;
    
    cont_timechunk=0; % counter of each time step of the chunk
    for i=IniLoop:length(tspan_chunk)  
      
      cont_timechunk=cont_timechunk+1;   
      
      % Save last position
      LastLL1 = LL1;
      LastLL2 = LL2;
      LastLL3 = LL3;      
      
      % Refloating
      if FateTypePart==1 & strcmpi(conf.Traj.Refloating,'yes')==1
        
        % Conditions for refloating: Monte Carlo Method and high tide (optional)
        if rand(1)<0.5^((tspan_chunk(i)-TimeLandPart)/conf.Traj.Tw)
          if strcmpi(conf.Traj.RefloatingAtHighTide,'no') || ~isempty(find(peak_pos==i))
            if conf.Traj.Verbose >= 2
              disp(['Refloating, particle: ' num2str(part) ', i:' num2str(i)])
            end
            % Refloating => Particle get its last known position in water
            LL1=PosWater1Part;
            LL2=PosWater2Part;
            LL3=PosWater3Part;
            LastLL1=PosWater1Part;
            LastLL2=PosWater2Part;
            LastLL3=PosWater3Part;
            FateTypePart = 0;
          end
        end
        
      % Resuspension
      elseif FateTypePart==2 & strcmpi(conf.Traj.Resuspension,'yes')==1
        
        Ush_partLonLat=interp3(longitude,latitude,TT,Ush,PosBottom1Part,PosBottom2Part,tspan_chunk(i));
        Vsh_partLonLat=interp3(longitude,latitude,TT,Vsh,PosBottom1Part,PosBottom2Part,tspan_chunk(i));
        
        % Transform into m/s
        [Ush_partLonLat,Vsh_partLonLat] = transformUScalar(PosBottom1Part,PosBottom2Part,Ush_partLonLat,Vsh_partLonLat,conf.OGCM.Coordinates,-1); 
        
        UshAbs = sqrt(Ush_partLonLat^2+Vsh_partLonLat^2);
        difDepth_partLonLat=interp3(longitude,latitude,TT,difDepth,PosBottom1Part,PosBottom2Part,tspan_chunk(i));
              
        if strcmpi(conf.Traj.KvOption,'fromOGCM')
          Kvsh_partLonLat=interp2(longitude,latitude,TT,Kvsh,PosBottom1Part,PosBottom2Part,tspan_chunk(i));
        else
          Kvsh_partLonLat=conf.Traj.Kv;
        end
        
        tau0=conf.Beh.WaterDensity*1000*Kvsh_partLonLat*UshAbs/difDepth_partLonLat;
               
        % Condition for resuspension
        if tau0 >= taucr2 
          if conf.Traj.Verbose >= 2
            disp(['Resuspension, particle: ' num2str(part) ', i:' num2str(i) ', Ush:' num2str(UshAbs) ', tau0:' num2str(tau0) ', taucr2:' num2str(taucr2)])
          end
          % Resuspension => Particle get its last known position 
          LL1=PosBottom1Part; 
          LL2=PosBottom2Part;
          LL3=PosBottom3Part; % Or option random value between bed and maximum height ?
          LastLL1=PosBottom1Part;
          LastLL2=PosBottom2Part;
          LastLL3=PosBottom3Part;
          FateTypePart = 0;
          
        % Condition for rolling/sliding
        elseif tau0 >= taucr1 & tau0 < taucr2 & strcmpi(conf.Traj.Sliding,'yes')
          % The acceleration can be taken into account for very small time steps (2DV for instance)
          % This is not appropriate for regional models, in this case, Up = conf.Beh.Cm*Ush
          %Up = LastUpPart + conf.Beh.Cm/conf.Beh.Cd*conf.Beh.WaterDensity/conf.Beh.ParticleDensity*conf.Beh.WaterViscosity/conf.Beh.ParticleSize*(TimeStepOut*86400)/difDepth_partLonLat*(Ush_partLonLat-LastUpPart);
          %Vp = LastVpPart + conf.Beh.Cm/conf.Beh.Cd*conf.Beh.WaterDensity/conf.Beh.ParticleDensity*conf.Beh.WaterViscosity/conf.Beh.ParticleSize*(TimeStepOut*86400)/difDepth_partLonLat*(Vsh_partLonLat-LastVpPart);
          Up = conf.Beh.Cm*Ush_partLonLat; 
          Vp = conf.Beh.Cm*Vsh_partLonLat;
          if conf.Traj.Verbose >= 2
            disp(['Rolling/sliding, particle: ' num2str(part) ', i:' num2str(i) ', Ush:' num2str(UshAbs) ', tau0:' num2str(tau0) ', taucr1:' num2str(taucr1)])
            disp(['particle: ' num2str(part) ', i:' num2str(i) ', Up: ' num2str(Up) ' Vp: ' num2str(Vp) ' LastUpPart: ' num2str(LastUpPart) ' LastVpPart: ' num2str(LastVpPart)])
          end         
          LastUpPart = Up;
          LastVpPart = Vp;
          
          % Transform into Lon/day and Lat/day before computing the movement
          [Up,Vp] = transformUScalar(PosBottom1Part,PosBottom2Part,Up,Vp,conf.OGCM.Coordinates,1);
          PosBottom1Part = PosBottom1Part + Up*TimeStepOut; 
          PosBottom2Part = PosBottom2Part + Vp*TimeStepOut; 
          
          if strcmpi(conf.OGCM.VerticalLayer,'sigma2depthVar') || strcmpi(conf.OGCM.BottomType,'Var') %Bottom varies with tide/time or morpho changes
            Depth_partLonLat=-interp3(longitude,latitude,TT,squeeze(Depth(:,:,end,:)),PosBottom1Part,PosBottom2Part,tspan_chunk(i));
          elseif strcmpi(conf.OGCM.VerticalLayer,'hybrid')
            Depth_partLonLat=interp3(longitude,latitude,TT,repmat(H,1,1,length(TT))+E,PosBottom1Part,PosBottom2Part,tspan_chunk(i));
          else
            Depth_partLonLat=griddata(Lon(mask_water==1),Lat(mask_water==1),H(mask_water==1),PosBottom1Part,PosBottom2Part,'nearest'); %#ok<GRIDD> %For the idealized bay
          end
          PosBottom3Part = -Depth_partLonLat;
          LL1=PosBottom1Part; 
          LL2=PosBottom2Part;
          LL3=PosBottom3Part;
          LastLL1=PosBottom1Part;
          LastLL2=PosBottom2Part;
          LastLL3=PosBottom3Part;
  
        % No motion  
        else
          if conf.Traj.Verbose >= 2
            disp(['No Resuspension, particle: ' num2str(part) ', i:' num2str(i) ', Ush:' num2str(UshAbs) ', tau0:' num2str(tau0) ', taucr1:' num2str(taucr1)])  
          end
          LastUpPart = 0.0;
          LastVpPart = 0.0;
        end
      end
        
      % Definition of internal time steps
      tspan=tspan_chunk(i)-direction*TimeStepOut:direction*TimeStepCalc:tspan_chunk(i);
      
      % Optimization of the interpolations by reducing the data matrices in the temporal dimension
      %   disp('tspan : ')
      %   disp(num2str(tspan))
      %   disp('TT : ')
      %   disp(num2str(TT))
      [~,Ind1] = min(abs(TT-tspan(1)));
      %   disp(['Ind1 : ' num2str(Ind1)])
      if direction > 0
        if TT(Ind1(1))>tspan(1)
          Ind1=Ind1-1;
        end
      else
        if TT(Ind1(1))<tspan(1)
          Ind1=Ind1-1;
        end  
      end
      %   disp(['Ind1 : ' num2str(Ind1)])
      [~,Ind2] = min(abs(TT-tspan(end)));
      %   disp(['Ind2 : ' num2str(Ind2)])
      if direction > 0
        if TT(Ind2(1))<tspan(end)
          Ind2=Ind2+1;
        end
      else
        if TT(Ind2(1))>tspan(end)
          Ind2=Ind2+1;
        end  
      end
      %   disp(['Ind2 : ' num2str(Ind2)])
      Ind = Ind1(1):Ind2(1);
      TTs = TT(Ind);  
      Us = U(:,:,:,Ind);
      Vs = V(:,:,:,Ind);
      Ws = W(:,:,:,Ind);
      Depths = Depth(:,:,:,Ind);
      KHs = 0.0;
      if strcmpi(conf.Traj.KhOption,'fromOGCM')
        KHs = KH(:,:,:,Ind);
      end
      KVs = 0.0;
      if strcmpi(conf.Traj.KvOption,'fromOGCM')
        KVs = KV(:,:,:,Ind);
      end
      
      % Particle in water => Advection and turbulence computation
      if FateTypePart==0
         
        for jj=2:length(tspan)
          
          % Advection
          if strcmpi(conf.Traj.Mode,'2D')
            [ln_ADV,lt_ADV] = AdvectionSchemes2D(conf.Traj.Scheme,longitude,latitude,TTs,Us,Vs,LL1,LL2,tspan(jj-1),TimeStepCalc,direction);
            h_ADV = 0.0;
          else
            [ln_ADV,lt_ADV,h_ADV] = AdvectionSchemes(conf.Traj.Scheme,longitude,latitude,Depths,TTs,Us,Vs,Ws,LL1,LL2,LL3,tspan(jj-1),TimeStepCalc,direction);
          end

          % CFL Check for advection, movement must be less than one grid cell
          CFLX = abs(ln_ADV-LL1)/DLon;
          CFLY = abs(lt_ADV-LL2)/DLat;
         if conf.Traj.Verbose >= 1
           if CFLX > 1 || CFLY > 1
             disp(['Warning, Advection CFL > 1 for part ' num2str(part) ', CFLX: ' num2str(CFLX,'%.2f') ', CFLY: ' num2str(CFLY,'%.2f')])
           end
         end
          
          % Turbulence
          if strcmpi(conf.Traj.KhOption,'fromOGCM')
            kh=interpDispersion(longitude,latitude,Depths,TTs,KHs,[LL1 LL2 LL3],tspan(jj)); %IJR April 2020
            % IJR22: To be checked!!!
          else
            kh=conf.Traj.Kh;
          end
          if strcmpi(conf.Traj.KvOption,'fromOGCM')
            kv=interpDispersion(longitude,latitude,Depths,TTs,KVs,[LL1 LL2 LL3],tspan(jj)); %IJR April 2020
          else
            kv=conf.Traj.Kv;
          end

          [TurbHx, TurbHy, TurbV] = HTurb30(1,TimeStepCalc*24*60*60,'Kh',kh,'Kv',kv);

          if ~strcmpi(conf.OGCM.Coordinates,'cartesian')
            dlnTur=km2deg(TurbHx/1000);
            dltTur=km2deg(TurbHy/1000);
          else
            dlnTur=TurbHx;
            dltTur=TurbHy;
          end

          % Behaviour/Sinking
          % Settling velocity at each time stamp
          if strcmpi(conf.Traj.Mode,'2D')
            TurbV = 0.0
            h_BEH = 0.0
          else
            Behaviour = behaviour(conf.Beh,tspan,ReleaseTime); %IJR 0424
			
			% Floating particles don't come from the bed even in backward computation
            if Behaviour.Ws(jj)>0 && direction==-1
                h_BEH = -Behaviour.Ws(jj)*(TimeStepCalc*24*60*60); % Ws (m/s) and TimeStepCalc (days) => h_BEH (m)
            else
                h_BEH = Behaviour.Ws(jj)*(TimeStepCalc*24*60*60); % Ws (m/s) and TimeStepCalc (days) => h_BEH (m)
            end
          end

          % New particule position
          LL1 = ln_ADV+dlnTur;     % Advection + Turbulence
          LL2 = lt_ADV+dltTur;     % Advection + Turbulence
          LL3 = h_ADV+TurbV-h_BEH; % Advection + Turbulence + Behaviour 
          
        end            
      end
      
      % 3D boundary condition check: 
      %   In water:      Fate==0
      %   Beaching:      Fate==1
      %   Deposition:    Fate==2
      %   Out of domain: Fate==3
                  
      % Domain check: 0 (out of domain) or 1 (inside of domain)
      inDomain=inpolygon(LL1,LL2,LonDomain,LatDomain);
      if inDomain==0
        if conf.Traj.Verbose >= 2
          disp(['Out of domain, particle: ' num2str(part) ', i:' num2str(i)])
        end
        LL1 = NaN;
        LL2 = NaN;
        LL3 = NaN;
        FateTypePart = 3;
        TimeOutDomainPart=tspan(end); %or tsapn_chunk(i) IJR 14/02/24
        
      % In domain
      else      
        inLand=inpolygon(LL1,LL2,xland,yland);
        
        % In land (but not already inland), beaching tests
        if inLand==1 & FateTypePart~=1
          if conf.Traj.Verbose >= 2
            disp(['InLand, particle: ' num2str(part) ', i:' num2str(i)])
          end
          TimeLandPart=tspan_chunk(i); %or tspan(end)
          
          % Avoid particles higher than water surface (even if in land)
          if LL3>0.0 % 0, 0.1 flotteurs
            LL3=0.0;
            if conf.Traj.Verbose >= 2
               disp(['Out of water, particle: ' num2str(part) ', i:' num2str(i)])
            end
          end

          % Save the last position of particles in the water (for refloating)
          if strcmpi(conf.Traj.Beaching,'yes')==1            
            FateTypePart = 1;
            Depth_partLonLat = 0.0;
            PosWater1Part=LastLL1; 
            PosWater2Part=LastLL2; 
            PosWater3Part=LastLL3;             
          else
            LL1 = LastLL1;
            LL2 = LastLL2;
            LL3 = LastLL3;
          end
        end 
              
        % Depth at particle position (to check deposition)
        if strcmpi(conf.OGCM.VerticalLayer,'sigma2depthVar') || strcmpi(conf.OGCM.BottomType,'Var') %Bottom varies with tide/time or morpho changes
          Depth_partLonLat=-interp3(longitude,latitude,TTs,squeeze(Depths(:,:,end,:)),LL1,LL2,tspan(end));
        elseif strcmpi(conf.OGCM.VerticalLayer,'hybrid')
          Depth_partLonLat=interp3(longitude,latitude,TT,repmat(H,1,1,length(TT))+E,LL1,LL2,tspan(end));
        else
          Depth_partLonLat=griddata(Lon(mask_water==1),Lat(mask_water==1),H(mask_water==1),LL1,LL2,'nearest'); %#ok<GRIDD> %For the idealized bay
        end
        
        % In water
        if FateTypePart==0
          
          % Avoid particles higher than water surface
          if LL3>0.0 % 0, 0.1 flotteurs
            LL3=0.0;
            if conf.Traj.Verbose >= 2
                disp(['Out of water, particle: ' num2str(part) ', i:' num2str(i)])
            end
          end
          
          % Deposition check
          if abs(LL3)>=Depth_partLonLat && strcmpi(conf.Traj.Mode,'3D')==1 % IJR 2D
            if conf.Traj.Verbose >= 2
              disp(['Reaching the bed, particle: ' num2str(part) ', i:' num2str(i) ', at depth:' num2str(Depth_partLonLat,'%.2f') ', LL3: ' num2str(LL3,'%.2f')])
            end
            if strcmpi(conf.Traj.Deposition,'yes')==1
              FateTypePart=2;
              PosBottom1Part=LL1;
              PosBottom2Part=LL2;
              PosBottom3Part=-Depth_partLonLat;        % => other option ? 
              LastUpPart = (LL1-LastLL1)/TimeStepOut;  % => 0 is another option
              LastVpPart = (LL2-LastLL2)/TimeStepOut;  % => 0 is another option
              [LastUpPart,LastVpPart] = transformUScalar(LL1,LL2,LastUpPart,LastVpPart,conf.OGCM.Coordinates,-1); % convert to m/s
              TimeSettlingPart=tspan(end);     %IJR 24/02/24        
            else
              LL3 = min([-Depth_partLonLat 0]); %min to avoid problems when interpolation near the coast
            end
          end
        end              
      end
      
      % Save the bottom depth at each particle position
      DepthBottomTrajPart(i)=Depth_partLonLat;     
      
      % Save trajectory position at output time step
      ln_Tot_Part(i)=LL1;
      lt_Tot_Part(i)=LL2;
      h_Tot_Part(i)=LL3;
      ts_Tot_Part(i)=tspan(end);
    end
    DepthBottomTraj(part,:) = DepthBottomTrajPart;
    ln_Tot_Output(part,:)=ln_Tot_Part;
    lt_Tot_Output(part,:)=lt_Tot_Part;
    h_Tot_Output(part,:)=h_Tot_Part;
    ts_Tot_Output(part,:)=ts_Tot_Part;
    FateType(part) = FateTypePart;  
    LastUpOutput(part) = LastUpPart;
    LastVpOutput(part) = LastVpPart;
    TimeLand(part) = TimeLandPart; %IJR 14/02/24
    TimeOutDomain(part) = TimeOutDomainPart; %IJR 14/02/24
    TimeSettling(part) = TimeSettlingPart; %IJR 14/02/24
    PosWater1(part) = PosWater1Part;
    PosWater2(part) = PosWater2Part;
    PosWater3(part) = PosWater3Part;
    PosBottom1(part) = PosBottom1Part;
    PosBottom2(part) = PosBottom2Part;
    PosBottom3(part) = PosBottom3Part;
    Duration(part) = toc;
  end
  
  % For parallel computing, the input arrays cannot be also the output arrays
  ln_Tot_Input=ln_Tot_Output(:,end);
  lt_Tot_Input=lt_Tot_Output(:,end);
  h_Tot_Input=h_Tot_Output(:,end);
  LastUpInput(:) = LastUpOutput(:);
  LastVpInput(:) = LastVpOutput(:);
  
  % Processing time informations
  if conf.Traj.Verbose >= 1
    disp(['Mean time to process each particle trajectory: ' num2str(mean(Duration),'%.1f') ' s'])
  end
     
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Save data and update data for the external loop
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % Build the structure TRAJ for each chunk
  [ln,lt,h] = deal(ln_Tot_Output,lt_Tot_Output,h_Tot_Output);
  TRAJ = TRAJstruct( size(ln) ); %call the strucutre TRAJ
  
  % Save particle positions
  TRAJ.Lon = ln;
  TRAJ.Lat = lt;
  TRAJ.Depth = h;
  
  TRAJ.TimeStamp = ts_Tot_Output(1,:); %Time %IJR22 (avant ts_Tot_Output(:)')
  TRAJ.InitialLonLatDepth = [TRAJ.Lon(:,1), TRAJ.Lat(:,1) TRAJ.Depth(:,1)]; %initial position
  TRAJ.FinalLonLatDepth(:,1) = TRAJ.Lon(:,end); %final position
  TRAJ.FinalLonLatDepth(:,2) = TRAJ.Lat(:,end);
  TRAJ.FinalLonLatDepth(:,3) = TRAJ.Depth(:,end);
  
  TRAJ.DepthBottomTraj = DepthBottomTraj; %Bottom depth along the particles trajectory
  TRAJ.FateType=FateType; %Fate type
  
  TRAJ.TimeSettling=TimeSettling; %IJR 14/02/24
  TRAJ.TimeOutDomain=TimeOutDomain; %IJR 14/02/24
  TRAJ.TimeBeaching=TimeLand; %IJR 14/02/24
  
  % Saving trajectory
  try
    if ~exist( conf.Traj.BaseDir, 'dir' )
      fprintf('%s is creating the directory:\n%s\n',mfilename,conf.Traj.BaseDir);
      mkdir(conf.Traj.BaseDir);
    end
    filename = fullfile(...
      conf.Traj.BaseDir,...
      strcat(conf.Traj.ScenarioName,'_Chunk',num2str(chunk),...
      '.mat'));
    save(filename,'TRAJ');
  catch
    fprintf('Saving trajectory failed because:\n')
    res = lasterror;
    fprintf('%s\n',res.message)
  end
  
  % Clear data
  clear TRAJ pom
  
  % Time per chunk
  ChunkDuration = 24*3600*(now-TheDate);
  disp(['chunk processing time: ' num2str(floor(ChunkDuration/60)) ' min ' num2str(mod(ChunkDuration,60),'%.1f') ' s'])
  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Join chunks and save results (Matlab format for the moment)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TRAJ = TRAJstruct([numpar length(tspanTOT)]) ;

% Join main data
inipos=1;
for i=1:chunk
  filename = fullfile(conf.Traj.BaseDir,strcat(conf.Traj.ScenarioName,...
    '_Chunk',num2str(i),'.mat'));
  TRAJ_chunk=load(filename);
  
  num_timestamps=length(TRAJ_chunk.TRAJ.TimeStamp);
  
  TRAJ.Lon(:,inipos:inipos+num_timestamps-1)=TRAJ_chunk.TRAJ.Lon;
  TRAJ.Lat(:,inipos:inipos+num_timestamps-1)=TRAJ_chunk.TRAJ.Lat;
  TRAJ.Depth(:,inipos:inipos+num_timestamps-1)=TRAJ_chunk.TRAJ.Depth;
  TRAJ.TimeStamp(:,inipos:inipos+num_timestamps-1)=TRAJ_chunk.TRAJ.TimeStamp;
  TRAJ.DepthBottomTraj(:,inipos:inipos+num_timestamps-1)=TRAJ_chunk.TRAJ.DepthBottomTraj;
  
  inipos=inipos+num_timestamps;
  
  delete(filename)
  
end

% Transform FateType code to character
TRAJ.FateType=cell(numpar,1);
TRAJ.FateType(TRAJ_chunk.TRAJ.FateType==0)={'Water',};
TRAJ.FateType(TRAJ_chunk.TRAJ.FateType==1)={'Land'};
TRAJ.FateType(TRAJ_chunk.TRAJ.FateType==2)={'Bottom'};
TRAJ.FateType(TRAJ_chunk.TRAJ.FateType==3)={'OutDomain'};

TRAJ.TimeSettling=TRAJ_chunk.TRAJ.TimeSettling;
TRAJ.TimeBeaching=TRAJ_chunk.TRAJ.TimeBeaching;
TRAJ.TimeOutDomain=TRAJ_chunk.TRAJ.TimeOutDomain;

clear TRAJ_chunk


% Initial and final positions and time of beaching, settling or outdomain
TRAJ.InitialLonLatDepth = [TRAJ.Lon(:,1), TRAJ.Lat(:,1) TRAJ.Depth(:,1)];

for j=1:numpar
  Finalpos = find(isnan(TRAJ.Lat(j,:))); %Last position in water
  aux=find(diff(Finalpos)~=1);
  
  %    if ~isempty(Finalpos) && Finalpos(1)~=1 && isempty(aux) % for particles going out of water 1 time
  if isnan(TRAJ.Lon(j,end)) && isempty(aux) && Finalpos(1)~=1
    TRAJ.FinalLonLatDepth(j,1) = TRAJ.Lon(j,Finalpos(1)-1);
    TRAJ.FinalLonLatDepth(j,2) = TRAJ.Lat(j,Finalpos(1)-1);
    TRAJ.FinalLonLatDepth(j,3) = TRAJ.Depth(j,Finalpos(1)-1);
    
    TRAJ.TrajectoryDuration(j)=abs(TRAJ.TimeStamp(Finalpos(1)-1)-TRAJ.TimeStamp(1)); %days
    
  elseif ~isempty(aux) && isnan(TRAJ.Lon(j,end)) % for particles replace in water several times
    
    TRAJ.FinalLonLatDepth(j,1) = TRAJ.Lon(j,Finalpos(aux(end)+1)-1);
    TRAJ.FinalLonLatDepth(j,2)= TRAJ.Lat(j,Finalpos(aux(end)+1)-1);
    TRAJ.FinalLonLatDepth(j,3)= TRAJ.Depth(j,Finalpos(aux(end)+1)-1);
    
    TRAJ.TrajectoryDuration(j)=abs(TRAJ.TimeStamp(Finalpos(aux(end)+1)-1)-TRAJ.TimeStamp(1)); %days
    
  else %particle in the water at the end of the simulation
    TRAJ.FinalLonLatDepth(j,1) = TRAJ.Lon(j,end);
    TRAJ.FinalLonLatDepth(j,2) = TRAJ.Lat(j,end);
    TRAJ.FinalLonLatDepth(j,3) = TRAJ.Depth(j,end);
    
    TRAJ.TrajectoryDuration(j)=abs(TRAJ.TimeStamp(end)-TRAJ.TimeStamp(1)); %days
    
  end
  
end

% Number of particles for each fate type
TRAJ.WaterParticles=length(find(strcmpi(TRAJ.FateType,'Water')));
TRAJ.BeachedParticles=length(find(strcmpi(TRAJ.FateType,'Land')));
TRAJ.BottomParticles=length(find(strcmpi(TRAJ.FateType,'Bottom')));
TRAJ.OutDomainParticles=length(find(strcmpi(TRAJ.FateType,'OutDomain')));

% Other metadata
TRAJ.conf = conf;
Behaviour = behaviour(conf.Beh,TRAJ.TimeStamp,ReleaseTime); %IJR Feb24
TRAJ.conf.Beh= Behaviour; %IJR Feb24

TRAJ.Domain.x=xland;
TRAJ.Domain.y=yland;
TRAJ.Grid=infogrid;


% Save
filename = fullfile(conf.Traj.BaseDir,strcat(conf.Traj.ScenarioName,'.mat'));
save(filename,'TRAJ');

fprintf('+++++++ Finished at: %s\n',datestr(now,0));
disp('++++++++++++ $$$$$$$$$$$$$$$$$$ Exit $$$$$$$$$$$$$$$$ ++++++++++++')
end
