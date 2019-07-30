function transformPOMinputs_2D(conf_name)

% TRANSFORM OGCM outputs (SARCCM POM version) TO TrackMPD FORMAT

level=1; %level at which velocities are taken

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (SARCCM POM version)
% POM output files
% hmin (depth for land)

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,w,E,depth,time,time_str

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: transformFVCOMoutputs July 2019 Z ijalonrojas $
%
% Copyright (C) 2017-2019 Isabel Jalon-Rojas and Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

conf=feval(conf_name);

POM_Prefix = conf.OGCM.POM_Prefix;
POM_Suffix = conf.OGCM.POM_Suffix;
fnames=getAllFiles(conf.OGCM.BaseDir,strcat(POM_Prefix,'*',POM_Suffix),true);

% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%% Define model parameters (SARCCM POM Version)

t0=conf.OGCM.t0; %POM reference time variable (first time step)
Hmin =conf.OGCM.Hmin; %H for land

TimeStamps = 1:numel(fnames);
nTimeStamps=length(TimeStamps);
tstep=conf.OGCM.TimeStep;


%% Read and save the Grid info    

ncfile = fnames{1};

Lat = double(ncread(ncfile,'latitude'));      % units = 'degrees_north'
Lon = double(ncread(ncfile,'longitude'));     % units = 'degrees_east'

lvl = double(ncread(ncfile,'level'));         % units = 'm'.
BottomDepth = double(ncread(ncfile,'depth'));           % units = 'm'

% % Sigma coordinate
% Z=[0;cumsum(lvl(1:end-1))];
% KB = size(Z,1);
% ZZ=-0.5*(Z(1:KB-1)+Z(2:KB));
% ZZ(KB)=2*ZZ(KB-1)-ZZ(KB-2);


% (water/land) mask
mask = double(BottomDepth~=1);
i=BottomDepth<Hmin;
BottomDepth(i)=Hmin;
mask2=~i;
mask=double(mask & mask2);
mask_water = mask;
mask_water(mask_water~=0)=1;
mask_land = ~mask_water;

% mask_land3D = mask_land;
% for i=1:length(lvl)-1
%     mask_land3D = cat(3,mask_land3D,mask_land);
% end

%BottomDepth(mask_land)=NaN;


%% Read and save the variables varing with time

[numlon,numlat] = size(BottomDepth);
%numlvl = length(lvl);

%zeros_matrix=zeros(numlon,numlat,1);

for n=1:nTimeStamps
    
    ncfile = fnames{TimeStamps(n)};
    
% Velocity and elevation
    
    u_aux = double(ncread(ncfile,'u-velocity'));      % units = 'm/s'
    u_aux=u_aux(:,:,level);
    v_aux = double(ncread(ncfile,'v-velocity'));      % units = 'm/s'
    v_aux=v_aux(:,:,level);
    E = double(ncread(ncfile,'elevation'));       % units = 'm'
    
    
%     %Include nans in land grid points
%     u_aux(mask_land3D) = NaN;
%     v_aux(mask_land3D) = NaN;
%     w(mask_land3D) = NaN;
%     E(mask_land) = NaN;
    
    % Conversion from Arakawa C-grid to A-grid (IJR 21/05/2018)
    u=cat(1,u_aux(2,:,:),(u_aux(2:end-1,:,:)+u_aux(3:end,:,:))/2,u_aux(end,:,:)); 
    v=cat(2,v_aux(:,2,:),(v_aux(:,2:end-1,:)+v_aux(:,3:end,:))/2,v_aux(:,end,:));

    % We can also use the functions, but I guess it will take more
    % computational time
%     u=pom_rho3u_3d(u_aux,0);
%     v=pom_rho3v_3d(v_aux,0);
        
%Time
    time = t0+(TimeStamps(n)-1)*tstep;
    time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
    timestamps(n)=time;
    
    fprintf('changing format for time %s\n',time_str);
    
    
    u=permute(u,[2,1]); %lat,lon
    v=permute(v,[2,1]); %lat,lon,
    E=permute(E,[2,1]);
    
    % From m/s to cm/s
    
    u=u*100;
    v=v*100;
   % w=w*100;
    
    save([conf.Data.BaseDir '\TrackMPDInput' num2str(n) '.mat'],'u','v','E','time','time_str');

    
end

BottomDepth=permute(BottomDepth,[2,1]); %lat,lon
mask_water=permute(mask_water,[2,1]); %lat,lon

%save grid and time stamps
save([conf.Data.BaseDir '\grid.mat'],'Lat','Lon','BottomDepth','mask_water');
save([conf.Data.BaseDir '\timestamps.mat'],'timestamps');
fprintf('saving grid and timestamps\n');

