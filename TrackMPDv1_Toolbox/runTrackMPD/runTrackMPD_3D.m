function a = runTrackMPD_3D(conf_name)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script run TrackMPD model (3D mode)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: runTrackMPD_3D.m 001 2018-02-12 10:55:10Z ef $
%        Las version: 2018-06-19 ijalonrojas
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas & Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('running TrackMPD...\n')

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
par(:,pX)=M(:,1); %Create variable "par" with particles information
par(:,pY)=M(:,2);
par(:,pZ)=M(:,3);

numpar=size(par,1); 
disp(['Number of particles = ' num2str(numpar) ' ; Release Time = ' datestr(ReleaseTime,'dd-mmm-yyyy HH:MM') ' ; ' conf.Traj.Direction ' trajectory']);
    
% Domain

domain=load(conf.Data.Domain);
xland=domain(:,1);
yland=domain(:,2);

% Hydrodynamic model grid

infogrid=load([conf.Data.BaseDir '\grid.mat']);

longitude=infogrid.Lon; latitude=infogrid.Lat; H=infogrid.BottomDepth;
[Lon,Lat]=meshgrid(longitude,latitude); mask_water=infogrid.mask_water;

        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define time span of the hydrodynamic and trajectory models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time span of the hydrodynamic inputs
TimeStepInp=conf.OGCM.TimeStep;
load([conf.Data.BaseDir '\timestamps.mat']); % load timestamps 
tspanhydro=timestamps; 

% time span of the trajctories
TimeStepOut=conf.Traj.TimeStep;
tspanTOT=t0:direction*TimeStepOut:tend;

% check if the trajectory period (spanTOT) is covered by tspanhydro
% IMPORTANT(***): The calculation of advection in the internal loop needs
% data from two aditional time steps (TimeStepOut) after tend

if direction==1 %forward
   
    if t0<TimeStepOut(1) || tend+TimeStepOut>tspanhydro(end) %(***) 1 more time steps
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

% Initialize particles position and fate
LL = par; %Initial position
FateType=zeros(numpar,1); %Initial region 0=water; 1=land; 2=bottom; 3=out of Domain


% Calculate time span for the whole trajectory period

 % Select good option values for tracking (advection) - otherwise output is junk.
abs_tol = 1.0e-3; % IJR: Checked through sensitiviy analysis
rel_tol = 1.0e-3; % IJR: Checked through sensitiviy analysis
maxstep = 0.1/24; % IJR: Checked through sensitiviy analysis
options = odeset('RelTol',rel_tol,'AbsTol',abs_tol,'MaxStep',maxstep);

% Settling velocity at each time stamp
Behaviour = behaviour_new(conf.beh,tspanTOT);
ws=Behaviour.Ws;

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


% Start extertal loop: chunk method
fprintf('Starting chunk method; number of partitions: %i\n',numpartition);
chunk=0; %Initialize chunk counter
cont=1; %Initialize counter for each time step

for k=1:numpartition % EXTERNAL LOOP
   
    chunk=chunk+1; %Chunk numeration follows the forward or backward direction
    
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
    
    OGCMFileNames=cell(NumFilesChunk,1);
    contFile=1;
    for i=posFile1:direction:posFileEnd
        OGCMFileNames(contFile)={[conf.Data.BaseDir '\TrackMPDInput' num2str(i) '.mat']};
        contFile=contFile+1;
    end   
    
    
    % Read OGCM history for the chunk k
    
    TT=[];
    U=[];
    V=[];
    W=[];
    Depth=[];
    E=[];
    for i=1:NumFilesChunk
        
        dataOGCM_aux=load(OGCMFileNames{i});
        TT=cat(1,TT,dataOGCM_aux.time);
        U=cat(4,U,dataOGCM_aux.u);
        V=cat(4,V,dataOGCM_aux.v);
        W=cat(4,W,dataOGCM_aux.w);
        
        if strcmpi(conf.OGCM.VerticalLayer,'sigma2depthVar')
            Depth=cat(4,Depth,dataOGCM_aux.depth);
        elseif i==1 && (strcmpi(conf.OGCM.VerticalLayer,'sigma2depthCte') || strcmpi(conf.OGCM.VerticalLayer,'Rectangular'))
            Depth=cat(4,Depth,dataOGCM_aux.depth);
        elseif strcmpi(conf.OGCM.VerticalLayer,'sigma2depthCte')==0 && strcmpi(conf.OGCM.VerticalLayer,'sigma2dephtVar')==0 && strcmpi(conf.OGCM.VerticalLayer,'Rectangular')==0
            error ('Invalid type of vertical layer')
        end
        
        if strcmpi(conf.Traj.Refloating,'yes')
            E=cat(3,E,dataOGCM_aux.E);
        end
        
    end
    
    
    fprintf('******* Current time: %s  === Processing release time: %s Chunk: %s \n', ...
            datestr(now,0), datestr(ReleaseTime,0), num2str(chunk));
     
    % Initialize Refloating variables
    if strcmpi(conf.Traj.Refloating,'yes') && strcmpi(conf.Traj.Beaching,'no')
        error('refloating cannot occurr withouth beaching')
   
    elseif strcmpi(conf.Traj.Refloating,'yes') && strcmpi(conf.Traj.Beaching,'yes')
            Elevation_ts=squeeze(E(conf.Traj.xcoord_elev,conf.Traj.xcoord_elev,:)); % Time series of elevation
            Elevation_TRAJ=interp1(TT,Elevation_ts,tspan_chunk); %Interpolate the elevation time series at external time
            peak_pos=findpeaks(Elevation_TRAJ); %Identify high tide
    end
    refloat=zeros(numpar,length(tspan_chunk)); %refloat matrix [particles, time] equal to 1 when a particles is refloated at a given time

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% INTERNAL LOOP (IJR 25/05/2018
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %Initialize variables 
    len_tspan_chunk=length(tspan_chunk);
    ln_TOT = nan(len_tspan_chunk,numpar); 
    lt_TOT = nan(len_tspan_chunk,numpar);
    h_TOT = nan(len_tspan_chunk,numpar);
    ts_TOT = nan(len_tspan_chunk,1);
    DepthBottomTraj = nan(numpar,1); 
    
    % First value = Release position
    if chunk==1
        ts_TOT(1)=ReleaseTime;
        ln_TOT(1,:)=par(:,1);
        lt_TOT(1,:)=par(:,2);
        h_TOT(1,:)=par(:,3);
        DepthBottomTraj(:,1)=griddata(Lon(mask_water==1),Lat(mask_water==1),H(mask_water==1),LL(:,1),LL(:,2),'nearest'); %Bottom Depth at the release point
        IniLoop=2;
    else
        IniLoop=1;
    end
  
    cont_timechunk=0; %counter of each time step of the chunk

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Internal loop: calculate position for each time step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i=IniLoop:length(tspan_chunk)
        
        cont=cont+1; %time counter for the whole duration of the tracking
        cont_timechunk=cont_timechunk+1;
        
        tspan=[tspan_chunk(i)-direction*TimeStepOut tspan_chunk(i) tspan_chunk(i)+1*direction*TimeStepOut]; %three time steps (error with two values, the result has a lenght higher than 2)
        
        %Refloating  
        if strcmpi(conf.Traj.Refloating,'yes') && ~isempty(find(peak_pos==i)) %if high tide
                         
               for r=1:numpar
                   
                   if isnan(LL(r,1)) && FateType(r)==1 && rand(1)<0.5^((tspan_chunk(i)-TimeLand(r))/conf.Traj.Tw) %if particle is beached and Monte Carlo Method
                                        
                       LL(r,1)=PosWater(r,1); %Refloating=Update the particle postion on water
                       LL(r,2)=PosWater(r,2);
                       LL(r,3)=PosWater(r,3);
                       refloat(r,i)=1;
                             
                   end
               end
        end
        
        
        % Advection
        
        %if strcmpi(conf.Traj.Windage,'no')

            [ln_ADV,lt_ADV,h_ADV,ts_a] = feval( 'particle_track_ode_grid_LonLat',Lon,Lat,Depth,...
                                   U,V,W,conf.OGCM.VerticalLayer,TT,tspan,LL,options);
        
%         elseif strcmpi(conf.Traj.Windage,'yes')
%         
% %             [ln_ADV,lt_ADV,h_ADV,ts_a] = feval( 'particle_track_ode_grid_LonLat_Wind',Lon,Lat,Depth,...
% %                                    U,V,W,Uwind,Vwind,DensityRate,SurfaceRate,TT,tspan,LL,options );
%         else 
%             error ('Windage should be equal to yes or no')
%         end
        
        % Turbulence
        [TurbHx, TurbHy, TurbV] = HTurb30(numpar,TimeStepOut*24*60*60,'Kh',conf.Traj.Kh,'Kv',conf.Traj.Kv);
                
        dlnTur=km2deg(TurbHx/1000);
        dltTur=km2deg(TurbHy/1000);
        
        % Behaviour/Sinking
        
        h_BEH=ws(cont)*(TimeStepOut*24*60*60); % m
        
        
        % TOTAL
        
        ln_TOT(i,:)= ln_ADV(2,:)+dlnTur'; %ADV+TUR+BEH
        lt_TOT(i,:)= lt_ADV(2,:)+dltTur'; %ADV+TUR+BEH
        h_TOT(i,:)= h_ADV(2,:)+TurbV'-direction*h_BEH; %ADV+TUR+BEH
        
        ts_TOT(i)= ts_a(2); % time
        
        
        % Boundary condition check: Identify beaching particles(Fate=1)
                                           % settled particles  (Fate=2)
                                           % particles out of domain (Fate=3)
                                           % other: particles in the water (Fate=0)
                                  % Avoid particles out of water because or turbulence/sinking
        for part=1:numpar
        
            % Domain check: 0 (out of domain) or 1 (inside of domain)
            inDomain=inpolygon(ln_TOT(i,part),lt_TOT(i,part),[min(longitude) max(longitude)],[min(latitude) max(latitude)]);
            
            % Inland check: 0 (water) or 1 (land)
            inLand=inpolygon(ln_TOT(i,part),lt_TOT(i,part),xland,yland);
            
            % Depth at particle position (to check deposition)
            %Depth_partLonLat=interp2(Lon,Lat,Depth,ln_TOT(i,part),lt_TOT(i,part)); %That's good to simulate beaches slope
            Depth_partLonLat=griddata( Lon(mask_water==1),Lat(mask_water==1),H(mask_water==1),ln_TOT(i,part),lt_TOT(i,part),'nearest'); %#ok<GRIDD> %For the idealized bay
        
            % If there is beaching and refloating: check beaching and save
            % refloating conditions
            if inLand==1 && refloat(part,i)==0 && strcmpi(conf.Traj.Beaching,'yes')
                ln_TOT(i,part)=NaN; lt_TOT(i,part)=NaN;  h_TOT(i,part)=NaN;
                FateType(part)=1;
                TimeLand(part)=ts_TOT(i);
                if chunk~=1 && cont_timechunk==1 %the last postion of the first time step of a chunk is equal to LL (IJR 26/06/2019)
                    PosWater(part,1)=LL(part,1); %Save the last postion of particles in the water, for refloating
                    PosWater(part,2)=LL(part,2); %Save the last postion of particles in the water, for refloating
                    PosWater(part,3)=LL(part,3); %Save the last postion of particles in the water, for refloating
                else
                    PosWater(part,1)=ln_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                    PosWater(part,2)=lt_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                    PosWater(part,3)=h_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                end
            
            % If no beaching: check if particles in land and put them in the water (last position in water)    
            elseif inLand==1 && refloat(part,i)==0 && strcmpi(conf.Traj.Beaching,'no')
                FateType(part)=0;
                if chunk~=1 && cont_timechunk==1 %the last postion of the first time step of a chunk is equal to LL (IJR 26/06/2019)
                    ln_TOT(i,part)=LL(part,1); %Save the last postion of particles in the water, for refloating
                    lt_TOT(i,part)=LL(part,2); %Save the last postion of particles in the water, for refloating
                    h_TOT(i,part)=LL(part,3); %Save the last postion of particles in the water, for refloating
                else
                    ln_TOT(i,part)=ln_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                    lt_TOT(i,part)=lt_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                    h_TOT(i,part)=h_TOT(i-1,part); %Save the last postion of particles in the water, for refloating
                end
            
            % Settling check    
            elseif abs(h_TOT(i,part))>=Depth_partLonLat
                ln_TOT(i,part)=NaN; lt_TOT(i,part)=NaN;  h_TOT(i,part)=NaN;
                FateType(part)=2;
            
            % Out of domain check
            elseif inDomain==0 && ~isnan(ln_TOT(i,part))
                ln_TOT(i,part)=NaN; lt_TOT(i,part)=NaN;  h_TOT(i,part)=NaN;
                FateType(part)=3;
                
            % Reflated particle
            elseif refloat(part,i)==1
                FateType(part)=0; %Water for reloating part
            end  
        
            
            % Avoid particles out of water because or turbulence/sinking
            if h_TOT(i,part)>=0
                
                h_TOT(i,part)=-0.05; %Minimum possible depth
            end
            
            % Save the bottom depth at each particle position
            DepthBottomTraj(part,i)=Depth_partLonLat; 
            
            
        end
        
        
        % Initialize  particles positions for the next internal step
        LL=[ln_TOT(i,:)' lt_TOT(i,:)' h_TOT(i,:)']; 
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Save data and update data for the external loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Build the structure TRAJ for each chunk
    [ln,lt,h] = deal(ln_TOT',lt_TOT',h_TOT');
    TRAJ = TRAJstruct( size(ln) ); %call the strucutre TRAJ
    
    % Save particle positions
    TRAJ.Lon = ln;
    TRAJ.Lat = lt;
    TRAJ.Depth = h;

    TRAJ.TimeStamp = ts_TOT(:)'; %Time
    TRAJ.InitialLonLatDepth = [TRAJ.Lon(:,1), TRAJ.Lat(:,1) TRAJ.Depth(:,1)]; %initial position
    TRAJ.FinalLonLatDepth(:,1) = TRAJ.Lon(:,end); %final position
    TRAJ.FinalLonLatDepth(:,2) = TRAJ.Lat(:,end);
    TRAJ.FinalLonLatDepth(:,3) = TRAJ.Depth(:,end);
    
    TRAJ.DepthBottomTraj = DepthBottomTraj; %Bottom depth along the particles trajectory
    TRAJ.FateType=FateType; %Fate type


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
        
    % Update data for the external loop
    LL=TRAJ.FinalLonLatDepth;
        
    % Clear data
    clear TRAJ pom

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

clear TRAJ_chunk


% Initial and final positions and time of beaching, settling or outdomain
TRAJ.InitialLonLatDepth = [TRAJ.Lon(:,1), TRAJ.Lat(:,1) TRAJ.Depth(:,1)];

TRAJ.TimeBeaching=nan(numpar,1);
TRAJ.TimeSettling=nan(numpar,1);
TRAJ.TimeOutDomain=nan(numpar,1);
for j=1:numpar 
    Finalpos = find(isnan(TRAJ.Lat(j,:))); %Last position in water
    aux=find(diff(Finalpos)~=1); 
    
%    if ~isempty(Finalpos) && Finalpos(1)~=1 && isempty(aux) % for particles going out of water 1 time
     if isnan(TRAJ.Lon(j,cont)) && isempty(aux) && Finalpos(1)~=1   
        TRAJ.FinalLonLatDepth(j,1) = TRAJ.Lon(j,Finalpos(1)-1); 
        TRAJ.FinalLonLatDepth(j,2) = TRAJ.Lat(j,Finalpos(1)-1);
        TRAJ.FinalLonLatDepth(j,3) = TRAJ.Depth(j,Finalpos(1)-1);
        
        if strcmpi(TRAJ.FateType(j),'Land')
            TRAJ.TimeBeaching(j)=TRAJ.TimeStamp(Finalpos(1)-1);
        elseif strcmpi(TRAJ.FateType(j),'Bottom')
            TRAJ.TimeSettling(j)=TRAJ.TimeStamp(Finalpos(1)-1);
        elseif strcmpi(TRAJ.FateType(j),'OutDomain')
            TRAJ.TimeOutDomain(j)=TRAJ.TimeStamp(Finalpos(1)-1);
        end
        
        TRAJ.TrajectoryDuration(j)=abs(TRAJ.TimeStamp(Finalpos(1)-1)-TRAJ.TimeStamp(1)); %days
        
     elseif ~isempty(aux) && isnan(TRAJ.Lon(j,cont)) % for particles replace in water several times
            
             TRAJ.FinalLonLatDepth(j,1) = TRAJ.Lon(j,Finalpos(aux(end)+1)-1); 
             TRAJ.FinalLonLatDepth(j,2)= TRAJ.Lat(j,Finalpos(aux(end)+1)-1);
             TRAJ.FinalLonLatDepth(j,3)= TRAJ.Depth(j,Finalpos(aux(end)+1)-1); 
             
            if strcmpi(TRAJ.FateType(j),'Land')
                TRAJ.TimeBeaching(j)=TRAJ.TimeStamp(Finalpos(aux(end)+1)-1);
            elseif strcmpi(TRAJ.FateType(j),'Bottom')
                TRAJ.TimeSettling(j)=TRAJ.TimeStamp(Finalpos(aux(end)+1)-1);
            elseif strcmpi(TRAJ.FateType(j),'OutDomain')
                TRAJ.TimeOutDomain(j)=TRAJ.TimeStamp(Finalpos(aux(end)+1)-1);
            end
            
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
TRAJ.behaviour= Behaviour;

TRAJ.Domain.x=xland;
TRAJ.Domain.y=yland;
TRAJ.Grid=infogrid;


% Save
filename = fullfile(conf.Traj.BaseDir,strcat(conf.Traj.ScenarioName,'.mat'));
save(filename,'TRAJ');

fprintf('+++++++ Finished at: %s\n',datestr(now,0));
disp('++++++++++++ $$$$$$$$$$$$$$$$$$ Exit $$$$$$$$$$$$$$$$ ++++++++++++')
end
