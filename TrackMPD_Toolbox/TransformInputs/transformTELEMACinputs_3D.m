function transformTELEMACinputs_3D(conf_name,confOGCM_name)
% TRANSFORM OGCM outputs (FURG TELEMAC version) to TrackMPD format
% Itele and I.Jalon-Rojas  12 Nov 2019; based on LoadFVCOMFiles_3D


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (TELEMAC version)
% TELEMAC output file
% numlat number of points in the new rectangular grid (latitude dimension)
% numlon number of points in the new rectangular gird (longitude dimension)
% 

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,w,E,depth,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

%conf=feval(conf_name); version 1
conf=feval(conf_name);
confOGCM=feval(confOGCM_name); %IJR new input format
conf=mergeStructure(conf,confOGCM); %IJR new input format

file = conf.OGCM.TELEMACFile;
domain=load(conf.Data.Domain);
t0=conf.OGCM.t0; %TELEMAC reference time variable (first time step)

% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%% Define model parameters 

numlon=conf.OGCM.NumLonGrid; 
numlat=conf.OGCM.NumLatGrid;


% Read TELEMAC data

data = telheadr(file); %le o resultado do telemac

NumGridPts = data.NPOIN/data.NPLAN; 
numlvl = data.NPLAN;
NTimeStamps = data.NSTEPS;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the Grid info
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = data.XYZ(1:NumGridPts,1); %1st layer (Superfície)
y = data.XYZ(1:NumGridPts,2); %1st layer (Superfície)  


%Define boundaries of the new grid
if strcmpi(conf.OGCM.cut,'yes')
    fprintf('cutting grid\n');
    minLon = conf.OGCM.minLon; % min lon
    maxLon = conf.OGCM.maxLon; % max lon
    minLat = conf.OGCM.minLat; % min lat
    maxLat = conf.OGCM.maxLat; % max lat
else
    minLon = min(x); % min lon
    maxLon = max(x); % max lon
    minLat = min(y); % min lat
    maxLat = max(y); % max lat
end


% Tranformation to rectangular grid

Lat=linspace(minLat,maxLat,numlat);
Lon=linspace(minLon,maxLon,numlon);

[Lon_matrix,Lat_matrix]=meshgrid(Lon,Lat);


ti=ones(size(Lon_matrix));

for i=1:numlat
    for j=1:numlon

        ti(i,j) = ~inpolygon(Lon(j),Lat(i),domain(:,1),domain(:,2)); %ti=0-->Land point

    end
end

ti(ti==0)=NaN;

% (water/land) mask
mask_water=zeros(size(ti));
mask_water(~isnan(ti))=1;
mask_land = ~mask_water;

mask_land3D = mask_land;
for i=1:numlvl-1
    mask_land3D = cat(3,mask_land3D,mask_land);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save Time information    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TT = 0:data.DT:data.NSTEPS*data.DT; %dt is in seconds
timestamps=t0+TT/60/60/24; %transform to days
save([conf.Data.BaseDir '/timestamps.mat'],'timestamps');
fprintf('saving timestamps\n');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the variables varing with time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read variables (NumGridPts,layer)
data_t=data;

for t = 1:NTimeStamps %carrega p/cada timestep
 
    % Open data for each time step
    data_t = telstepr(data_t,t); %carrega p/cada timestep   
    
    %depth unit: m
    for lv=1:numlvl
        DEPTH(:,lv,t) = data_t.RESULT(NumGridPts*(numlvl-lv)+1:NumGridPts*(numlvl+1-lv),1);
        UU(:,lv,t) = data_t.RESULT(NumGridPts*(numlvl-lv)+1:NumGridPts*(numlvl+1-lv),2);
        VV(:,lv,t) = data_t.RESULT(NumGridPts*(numlvl-lv)+1:NumGridPts*(numlvl+1-lv),3);
        %velocities units = 'm/s'
    end
    
end

% Calculate Bottom Depth from Depth (last layer, first time step)

BottomDepth=griddata(x,y,DEPTH(:,numlvl,1),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of Depth in the new grid
BottomDepth(mask_water==0)=NaN; %Land point=0
%BottomDepth(isnan(BottomDepth))=1; 
%BottomDepth=permute(BottomDepth,[2,1]); %lat lon

% save grid
save([conf.Data.BaseDir '/grid.mat'],'Lat','Lon','mask_water','BottomDepth');
fprintf('saving grid\n');

% Loop for each time step

%zeros_matrix=zeros(numlat,numlon,1);

for i=1:NTimeStamps
    
    depth=nan(numlat,numlon,numlvl);
    u=nan(numlat,numlon,numlvl);
    v=nan(numlat,numlon,numlvl);
    w=nan(numlat,numlon,numlvl);
    E=nan(numlat,numlon);
    BottomDepth=nan(numlat,numlon);
    
    % 3D variables
    for j=1:numlvl
        Uaux=griddata(x,y,UU(:,j,i),Lon_matrix,Lat_matrix,'linear'); % Interpolation of U in the new grid
        Uaux(mask_water==0)=0; %Land point=0
        Uaux(isnan(Uaux))=0; 
        u(:,:,j)=Uaux;

        Vaux=griddata(x,y,VV(:,j,i),Lon_matrix,Lat_matrix,'linear'); % Interpolation of U in the new grid
        Vaux(mask_water==0)=0;
        Vaux(isnan(Vaux))=0; 
        v(:,:,j)=Vaux;
        
        Depth_aux=griddata(x,y,DEPTH(:,j,i),Lon_matrix,Lat_matrix,'linear'); % Interpolation of U in the new grid
        Depth_aux(mask_water==0)=NaN; %Land point=NaN
       % Depth_aux(isnan(Depth_aux))=1; 
        depth(:,:,j)=Depth_aux;
        
        clear Uaux Vaux Depth_aux
    end
    
    w=zeros(size(u));
    
    %2D variables
    
    E=depth(:,:,1);
    BottomDepth=depth(:,:,end);
    
    
    %1D variables: time
    time=timestamps(i);
    time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
    fprintf('changing format for time %s\n',time_str);
    
    
    % From m/s to cm/s 
    u=u*100;
    v=v*100;
    w=w*100;
    
    %Change in the reference system (surface constant, varying bottom)
    depth=depth-repmat(depth(:,:,1),[1,1,size(depth,3)]); %Change in the reference system (surface constant, varying bottom)
    
    % Land point ==1 for TrackMPD
    depth(isnan(depth))=1; 
    
    % change dimensions from lon,lat,z to lat,lon,z
%     u=permute(u,[2,1,3]); %lat,lon,z
%     v=permute(v,[2,1,3]); %lat,lon,z
%     w=permute(w,[2,1,3]); %lat,lon,z
%     depth=permute(depth,[2,1,3]); %lat,lon
%     E=permute(E,[2,1]);
    
    for ii=1:numlat
        for jj=1:numlon
            if sum(depth(ii,jj,:))==0
                depth(ii,jj,:)=1;
            end
        end
    end



    % save data for each time step 
    save([conf.Data.BaseDir '/TrackMPDInput' num2str(i) '.mat'],'u','v','w','E','time','time_str','depth','BottomDepth');
    
end

end
