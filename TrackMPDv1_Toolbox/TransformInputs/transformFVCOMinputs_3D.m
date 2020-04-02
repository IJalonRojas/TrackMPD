function transformFVCOMoutputs_3D(conf_name)
% TRANSFORM OGCM outputs (SARCCM FVCOM version) to TrackMPD format
% I.Jalon-Rojas  8 July 2019; based on LoadFVCOMFiles_3D


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (SARCCM FVCOM version)
% FVCOM output file
% numlat number of points in the new rectangular grid (latitude dimension)
% numlon number of points in the new rectangular gird (longitude dimension)
% name of time variable

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,w,E,depth,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Call the model configuration and inputs files

conf=feval(conf_name);

file = conf.OGCM.FVCOMFile;
grid = conf.OGCM.FVCOMGrid;

% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)


%% Define model parameters (SARCCM FVCOM Version)

numlon=conf.OGCM.NumLonGrid; 
numlat=conf.OGCM.NumLatGrid;
time_name=conf.OGCM.NameTime;
w_name=conf.OGCM.NameW;

[NumGridPts,numlvl,NTimeStamps]=size(ncread(file,'u'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the Grid info from Casename_Geo.grd file ZC  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% --load   *geo.grd file for grid information-----------
fid=fopen(grid);
[dims]=textscan(fid,'%d%d',1,'Headerlines',1);%Read header
[info_node]=textscan(fid,'%d%f%f%f',dims{1,2}(1)); %Read info of nodes
[info_cell]=textscan(fid,'%d%d%d%d%d',dims{1,1}(1)); %Read info of cells
cell_num=dims{1,1}(1);
lat_node=info_node{1,3};
lon_node=info_node{1,2};
cell_3nodeid(:,1)=info_cell{1,3};
cell_3nodeid(:,2)=info_cell{1,4};
cell_3nodeid(:,3)=info_cell{1,5};
%----calculate latitude and longtitude for each cell----
for idt_cell=1:cell_num
lat_cell(idt_cell)=lat_node(cell_3nodeid(idt_cell,1))+lat_node(cell_3nodeid(idt_cell,2))+lat_node(cell_3nodeid(idt_cell,3));
lat_cell(idt_cell)=lat_cell(idt_cell)/3;
lon_cell(idt_cell)=lon_node(cell_3nodeid(idt_cell,1))+lon_node(cell_3nodeid(idt_cell,2))+lon_node(cell_3nodeid(idt_cell,3));
lon_cell(idt_cell)=lon_cell(idt_cell)/3;
end

lat_v=lat_cell; % lat at cells
lon_v=lon_cell; %lon at cells
x=lon_node; %lon at nodes
y=lat_node; %lat at nodes

% For FVCOM using spherical coordinates

% lat_v=double(ncread(file,'latc')); % lat at cells
% lon_v=double(ncread(file,'lonc')); %lon at cells
% x=double(ncread(file,'lon')); %lon at nodes
% y=double(ncread(file,'lat')); %lat at nodes

nv=double(ncread(file,'nv')); %number of cells
TR = triangulation(nv,x,y); %triangular grid

% Tranformation to rectangular grid

Lat=linspace(min(lat_v),max(lat_v),numlat);
Lon=linspace(min(lon_v),max(lon_v),numlon);

[Lon_matrix,Lat_matrix]=meshgrid(Lon,Lat);

for i=1:numlat
    for j=1:numlon

        ti(i,j) = pointLocationQuadTree(TR,[Lon_matrix(i,j),Lat_matrix(i,j)]); %ti=Nan-->Land point

    end
end

% (water/land) mask
mask_water=zeros(size(ti));
mask_water(~isnan(ti))=1;
mask_land = ~mask_water;

mask_land3D = mask_land;
for i=1:numlvl-1
    mask_land3D = cat(3,mask_land3D,mask_land);
end

%Bottom Depth
h=double(ncread(file,'h'));
BottomDepth=double(griddata(x,y,h,Lon_matrix,Lat_matrix,'nearest')); % Interpolation of U in the new grid
BottomDepth(mask_water==0)=1; %Land=1

% Verification plot
figure;
mesh(Lon,Lat,BottomDepth);
hold on
plot3(x,y,h,'o');
title('Verification plot for grid transformation: Bottom depth in the new grid')

% save grid
save([conf.Data.BaseDir '\grid.mat'],'Lat','Lon','BottomDepth','mask_water');
fprintf('saving grid\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save Time information    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TT=double(ncread(file,time_name)); %modified julian date
timestamps=TT+678942; %all matlab versions
%TT=datetime(TT(:),'convertfrom','modifiedjuliandate'); % recent matlab version
%timestamps=datenum(TT); % recent matlab version
save([conf.Data.BaseDir '\timestamps.mat'],'timestamps');
fprintf('saving timestamps\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save the variables varing with time
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Read variables

UU=double(ncread(file,'u')); % units = 'm/s'
VV=double(ncread(file,'v')); % units = 'm/s'
WW=double(ncread(file,w_name)); % variable name of the vertical velocity from FVCOM

ele=double(ncread(file,'zeta')); %elevation (m)

siglay=double(ncread(file,'siglay')); 
h_con=double(griddata(x,y,h,Lon_matrix,Lat_matrix,'nearest')); %


% Loop for each time step

zeros_matrix=zeros(numlat,numlon,1);
depth=nan(numlat,numlon,numlvl);

for i=1:NTimeStamps
    
    u=nan(numlat,numlon,numlvl);
    v=nan(numlat,numlon,numlvl);
    w=nan(numlat,numlon,numlvl);
    E=nan(numlat,numlon);
    
    % 3D variables
    for j=1:numlvl
        Uaux=griddata(lon_v,lat_v,UU(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Uaux(mask_water==0)=0; %Land point=0
        Uaux(isnan(Uaux))=0; 
        u(:,:,j)=Uaux;

        Vaux=griddata(lon_v,lat_v,VV(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Vaux(mask_water==0)=0;
        Vaux(isnan(Vaux))=0; 
        v(:,:,j)=Vaux;
        
        Waux=griddata(lon_v,lat_v,WW(:,j,i),Lon_matrix,Lat_matrix,'nearest'); % Interpolation of U in the new grid
        Waux(mask_water==0)=0; 
        Waux(isnan(Waux))=0; 
        w(:,:,j)=Waux;
        
        %if i==1 %depth does not change with time in FVCOM (independent of elevation)
%         Depth_aux=siglayzc(:,:,j);
%         Depth_aux(mask_water==0)=0; 
%         depth(:,:,j)=Depth_aux;
%         end
        
        clear Uaux Vaux Waux
    end
    
    %2D variables: elevation
    ELEaux=double(griddata(x,y,ele(:,i),Lon_matrix,Lat_matrix,'nearest')); % Interpolation of U in the new grid
    ELEaux(mask_water==0)=0; %Land point=0
    ELEaux(isnan(ELEaux))=0; 
    E(:,:)=ELEaux;
   
    
    clear ELEaux
      
    %1D variables: time
    time=timestamps(i);
    time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
    fprintf('changing format for time %s\n',time_str);
    
    
    % From m/s to cm/s
    
    u=u*100;
    v=v*100;
    w=w*100;
    
    
    % Depth calculation
    sigma=[0 siglay(1,:) -1];
    BottomDepth(BottomDepth==0)=NaN;
      for sig=1:numlvl+2
        depth(:,:,sig)=(BottomDepth+E).*sigma(1,sig); %+ele-ele (-ele to make depth=0 at surface waters)
      end
      
      depth(isnan(depth))=1;
    
    
    u=cat(3,u(:,:,1),u,zeros_matrix);
    v=cat(3,v(:,:,1),v,zeros_matrix);
    w=cat(3,w(:,:,1),w,zeros_matrix);
    
    % save data for each time step 
    save([conf.Data.BaseDir '\TrackMPDInput' num2str(i) '.mat'],'u','v','w','E','time','time_str','depth');
    
end

end
