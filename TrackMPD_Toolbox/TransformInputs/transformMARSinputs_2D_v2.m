function transformMARSinputs_2D_v2(conf_name,confOGCM_name)
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

fileU = conf.OGCM.MARSFileU;
fileV = conf.OGCM.MARSFileV;
domain=load(conf.Data.Domain);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the Grid info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Arakawa-C grid
lat = ncread(fileU,'latitude');
lon = ncread(fileU,'longitude');

lat_u = ncread(fileU,'latitude_u');
lat_v = ncread(fileV,'latitude_v');

lon_u = ncread(fileU,'longitude_u');
lon_v = ncread(fileV,'longitude_v');

% dif_lat=min(diff(lat(1,:)));
% dif_lon=min(diff(lon(:,1)));

Lat=lat(1,:);
Lon=lon(:,1)';

%Define boundaries of the new grid
if strcmpi(conf.OGCM.cut,'yes')
    minLon = conf.OGCM.minLon; % min lon
    maxLon = conf.OGCM.maxLon; % max lon
    minLat = conf.OGCM.minLat; % min lat
    maxLat = conf.OGCM.maxLat; % max lat

    Lat=Lat(Lat>=minLat & Lat<=maxLat);
    Lon=Lon(Lon>=minLon & Lat<=maxLon);
    
end


BottomDepth = ncread(fileU,'H0')'; % units = 'm'
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

timestamps=ncread(fileU,'time')/(3600*24)+datenum(1900,1,1,0,0,0);
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
    
    u = ncread(fileU,'UZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    v = ncread(fileV,'VZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    %w = ncread(file,'WZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
    w = zeros(size(u));      % units = 'm/s'  
    E = ncread(fileU,'XE',[1 1 n],[Inf Inf 1]);       % units = 'm'    
     

% transform velocities
  
    % a) take surface velocities
    
    u=squeeze(u(:,:,end));    v=squeeze(v(:,:,end)); 
    
    % b) velocity at land = 0 m/s
    
     u(isnan(u))=0;
     v(isnan(v))=0;
    
    % c) interpolation at lat lon
    
     u=griddata(lat_u(1,:),lon_u(:,310),u,lat(1,:),lon(:,1)); 
     v=griddata(lat_v(1,:),lon_v(:,310),v,lat(1,:),lon(:,1));
    %!!!!!!!!!!!!!!!!!!!! % 80 only for ARCACHON application !!!!!!!!!!!!!!!!!
    
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
