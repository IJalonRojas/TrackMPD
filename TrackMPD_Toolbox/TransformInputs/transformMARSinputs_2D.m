function transformMARSinputs_2D(conf_name,confOGCM_name)
% TRANSFORM OGCM outputs (MARS-3D) to TrackMPD format (2D mode)
% I.Jalon-Rojas  20 Mars 2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (MARS version)
% MARS output file
% 

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,E,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

conf=feval(conf_name);
confOGCM=feval(confOGCM_name); %IJR new input format
conf=mergeStructure(conf,confOGCM); %IJR new input format

file = conf.OGCM.MARSFile;
domain=load(conf.Data.Domain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the Grid info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arakawa-C grid
lat = ncread(file,'latitude');
lon = ncread(file,'longitude');

lat_u = ncread(file,'latitude_u');
lat_v = ncread(file,'latitude_v');

lon_u = ncread(file,'longitude_u');
lon_v = ncread(file,'longitude_v');

% dif_lat=min(diff(lat(1,:)));
% dif_lon=min(diff(lon(:,1)));

Lat=lat(1,:);
Lon=lon(:,1)';

%Define boundaries of the new grid
if strcmpi(conf.OGCM.cut,'yes')
    
    fprintf('cutting grid\n');
    minLon = conf.OGCM.minLon; % min lon
    maxLon = conf.OGCM.maxLon; % max lon
    minLat = conf.OGCM.minLat; % min lat
    maxLat = conf.OGCM.maxLat; % max lat

    posLat=find(Lat>=minLat & Lat<=maxLat);
    posLon=find(Lon>=minLon & Lon<=maxLon);
    
    Lat=Lat(posLat);
    Lon=Lon(posLon);
    
    pos1Lat=posLat(1); posEndLat=posLat(end);
    pos1Lon=posLon(1); posEndLon=posLon(end);
else
    pos1Lat=1; posEndLat=length(Lat);
    pos1Lon=1; posEndLon=length(Lon);
    
end


BottomDepth = ncread(file,'H0')'; % units = 'm'
BottomDepth = BottomDepth(pos1Lat:posEndLat,pos1Lon:posEndLon);
%BottomDepth(isnan(BottomDepth))=1;

% (water/land) mask
mask_water=zeros(size(BottomDepth));
mask_water(~isnan(BottomDepth))=1;
mask_land = ~mask_water;

% save grid
save(fullfile(conf.Data.BaseDir,'grid.mat'),'Lat','Lon','BottomDepth','mask_water');
fprintf('saving grid\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save Time information    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

timestamps=ncread(file,'time')/(3600*24)+datenum(1900,1,1,0,0,0);
save(fullfile(conf.Data.BaseDir,'timestamps.mat'),'timestamps');
fprintf('saving timestamps\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the variables varing with time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% dimensions
[numlon,numlat] = size(BottomDepth);
nTimeStamps=length(timestamps);


for n=1:nTimeStamps %loop for each time step

% time    
    time = timestamps(n);
    time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
    fprintf('changing format for time %s\n',time_str);
    
% read velocity and elevation at each time step
    
    u = ncread(file,'UZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    v = ncread(file,'VZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    %w = ncread(file,'WZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    w = zeros(size(u));      % units = 'm/s'  
    E = ncread(file,'XE',[1 1 n],[Inf Inf 1]);       % units = 'm'    
     

% transform velocities
  
    % a) take surface velocities
    
    u=squeeze(u(:,:,end));    v=squeeze(v(:,:,end)); 
    
    % b) velocity at land = 0 m/s
    
     u(isnan(u))=0;
     v(isnan(v))=0;
    
    % c) interpolation at lat lon
    
     u=griddata(lat_u(1,:),lon_u(:,310),u,Lat',Lon); 
     v=griddata(lat_v(1,:),lon_v(:,310),v,Lat',Lon);
    %!!!!!!!!!!!!!!!!!!!! % 80 only for ARCACHON application !!!!!!!!!!!!!!!!!
     
    E=E(pos1Lon:posEndLon,pos1Lat:posEndLat);
    
    % d) permute: dim=latxlonxz
    
    u=permute(u,[2,1]); %lat,lon dimensions
    v=permute(v,[2,1]); %lat,lon dimensions
    E=permute(E,[2,1]); %lat, lon dimensions
     
    
    % d) from m/s to cm/s
    
     u=u*100;
     v=v*100;
    
%     % Dispersion coefficient 4/6/20
%     if strcmpi(conf.Traj.KhOption,'fromOGCM')
%            kv = ncread(file,'KZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
%            kv=squeeze(kv(:,:,end));
%            kv(isnan(kv))=0;
%            kv=permute(kv,[2,1]); %lat,lon dimensions
%         
%     end
     
    save(fullfile(conf.Data.BaseDir,['TrackMPDInput' num2str(n) '.mat']),'u','v','E','time','time_str');
    
end
