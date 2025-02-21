function transformMARSinputs_3D(conf_name,confOGCM_name)
% TRANSFORM OGCM outputs (MARS-3D) to TrackMPD format
% I.Jalon-Rojas  20 Mars 2020


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (MARS version)
% MARS output file
%

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,w,E,depth,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

conf=feval(conf_name);
confOGCM=feval(confOGCM_name); %IJR new input format
conf=mergeStructure(conf,confOGCM); %IJR new input format

file = conf.OGCM.MARSFile;
domain=load(conf.Data.Domain);

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

lvl = ncread(file,'level');   %sigma layer: level=SIG
%sig = ncread(file,'SIG');
lvl=flip(lvl); % flip vector: surface at line 1 and bottom at the last line

BottomDepth = ncread(file,'H0')'; % units = 'm'
BottomDepth = BottomDepth(pos1Lat:posEndLat,pos1Lon:posEndLon);

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
numlvl = length(lvl); nTimeStamps=length(timestamps);

ZZ=[0;lvl;-1]; %Add a layer for surface z=0 (depth(:,:,1)==E) and bottom
%to avoid problem when interpolating velocities near the surface and bottom


for n=1:nTimeStamps %loop for each time step

  % time
  time = timestamps(n);
  time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
  fprintf('changing format for time %s\n',time_str);

  % read velocity and elevation at each time step

  U = ncread(file,'UZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
  V = ncread(file,'VZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
  %w = ncread(file,'WZ',[1 1 1 n],[Inf Inf Inf 1]);      % units = 'm/s'
  W = zeros(size(U));      % units = 'm/s'
  E = ncread(file,'XE',[1 1 n],[Inf Inf 1]);       % units = 'm'


  % transform velocities

  % a) flip matrix: surface velocities at first line, bottom velocites at
  % the last line

  U=flip(U,3); V=flip(V,3); W=flip(W,3);

  % b) velocity at land = 0 m/s

  U(isnan(U))=0;
  V(isnan(V))=0;
  W(isnan(W))=0;

  % c) Interpolate at lon,lat

  for lv=1:numlvl
    u(:,:,lv)=griddata(lat_u(1,:),lon_u(:,310),U(:,:,lv),Lat',Lon);
    v(:,:,lv)=griddata(lat_v(1,:),lon_v(:,310),V(:,:,lv),Lat',Lon);
  end

  w=W(pos1Lon:posEndLon,pos1Lat:posEndLat,:);
  E=E(pos1Lon:posEndLon,pos1Lat:posEndLat);

  % d) permute: dim=latxlonxz

  u=permute(u,[2,1,3]); %lat,lon,z
  v=permute(v,[2,1,3]); %lat,lon,z
  w=permute(w,[2,1,3]); %lat,lon,z

  % e) Velocity at the new first layer (depth=0) = velocity at the old
  % first sigma layer
  % Velocity at the bottom = 0

  % Create a zero matrix with dimensions (lonxlat)
  zeros_matrix=zeros(numlon,numlat,1);

  u=cat(3,u(:,:,1),u,zeros_matrix);
  v=cat(3,v(:,:,1),v,zeros_matrix);
  w=cat(3,w(:,:,1),w,zeros_matrix);
  E=permute(E,[2,1]);

  %f) from m/s to cm/s

  u=u*100;
  v=v*100;
  w=w*100;

  % create depths at each grid point

  depth=nan(numlon,numlat,numlvl+2);
  for i=1:length(ZZ)
    % OLD REFERENCE SYSTEM REMOVED BY V. MARIEU, 2024/09/10
    depth_aux = ZZ(i)*(BottomDepth+E)+E; %+E-E; Tranform: sigma*(BottomDepth+E)+E
    depth_aux(isnan(depth_aux))= griddata(x,y,DEPTH(:,j,i),Lon_matrix(isnan(depth_aux)),Lat_matrix(isnan(depth_aux)),'nearest');
    depth(:,:,i)=depth_aux;
    % depth varying between -Bottomdepth and +E
    % -E: This is to change the reference system:
    %Surface: depth=0, bottom changing with tide
  end

  % MODIFIED BY MARIEU 2025/01 for interpolation close to the shore
  %depth(isnan(depth))=1;
  %depth(isnan(depth))=1
  %depth(depth>-0.01 & depth<0)=0; % To avoid problem with shallow waters (TEST!!!)


  save(fullfile(conf.Data.BaseDir,['TrackMPDInput' num2str(n) '.mat']),'u','v','w','E','time','time_str','depth');

  clear u v w E
end
