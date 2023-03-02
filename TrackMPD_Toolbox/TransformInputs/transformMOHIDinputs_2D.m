function transformMOHIDinputs_2D(conf_name,confOGCM_name) 
% TRANSFORM OGCM outputs (SARCCM FVCOM version) to TrackMPD format
% M. Desai and I.Jalon-Rojas 16 Sept 2019 modification of I.Jalon-Rojas  8 July 2019; based on LoadFVCOMFiles_3D
% 2D version by I.Jalon Rojas April 2020

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUTS (MOHID version)
% MOHID output files
% SigmaFile which has sigma layer stuff
% t0: first time step
% steps: number of time steps per MOHID output file

%%%% OUTPUTS
% grid.mat
% timestamps.mat
% One file for each time step containing u,v,E,time,time_str
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

conf=feval(conf_name);
confOGCM=feval(confOGCM_name); %IJR new input format
conf=mergeStructure(conf,confOGCM); %IJR new input format

MOHID_Prefix = conf.OGCM.MOHID_Prefix;
MOHID_Suffix = conf.OGCM.MOHID_Suffix;
fnames=getAllFiles(conf.OGCM.BaseDir,strcat(MOHID_Prefix,'*',MOHID_Suffix),true);

% sigma_file_path = conf.OGCM.SigmaFile;
% sigma_file =  fopen(sigma_file_path,'r');
% formatSpec = '%f';
% sigma = double(fscanf(sigma_file,formatSpec));
% ZZ = cumsum(sigma);
% ZZ=[0;ZZ]; %Add a layer for surface z=0 (depth(:,:,1)==E) to avoid 
%    %problem when interpolating velocities near the surface

    	

% TO AVOID MEMORY PROBLEMS: we will save one output file with the new format 
% for each OGCM model time step (the OGCM model time step is defined inside conf)

t0=conf.OGCM.t0; 

%Dealing with files:
nFiles = length(fnames);

%% Read and save the Grid info
ncfile = fnames{1};



fLat = double(ncread(ncfile,'lat'));      % units = 'degrees_north'
fLon = double(ncread(ncfile,'lon'));     % units = 'degrees_east'

%lvl = double(ncread(ncfile,'depth'));

%x = double(ncread(ncfile,'line_c'));
%y = double(ncread(ncfile,'column_c'));

%Lat = fLat; % IJR 25/05/20 mateus conf
%Lon = fLon; % IJR 25/05/20 mateus conf
Lat = fLat(1,:);
Lon = fLon(:,1).';

%%% Define boundaries of the new grid
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
%%%

numlat = length(Lat);
numlon = length(Lon);
%numlvl = length(lvl);
%zeros_matrix=zeros(numlat,numlon,1);

% % create constant bottom 
% num_cartesian_layers = 18;
% cartesian_depth =double(zeros(numlat,numlon,num_cartesian_layers));
% for j=1:num_cartesian_layers
% 	depth_lookup = num_cartesian_layers+1 -j; % depth variable has deeper depths indexed first
%    	for x_it=1:numlat
%    		for y_it = 1:numlon
%    			cartesian_depth(x_it,y_it,j) = lvl(depth_lookup);
%    		end
%    	end
% end

%%% note for tomrorow: some problem here


% (water/land) mask 

mask_water = double(ncread(ncfile,'mask'));
mask_water = mask_water(:,:,end).'; % surface
%mesh(fLon.',fLat.',mask_water); % show the watermask for sanity check '
mask_water = mask_water(pos1Lat:posEndLat,pos1Lon:posEndLon);
mask_land = ~mask_water;

BottomDepth = double(ncread(ncfile,'bathymetry'));
BottomDepth = BottomDepth.'; %'
BottomDepth = BottomDepth(pos1Lat:posEndLat,pos1Lon:posEndLon);
BottomDepth(mask_water==0)=NaN; %%IJR 030420


% Verification plot
%figure;
%mesh(fLon.',fLat.',BottomDepth); 
%hold on
%plot3(x,y,BottomDepth,'o');
%title('Verification plot for grid transformation: Bottom depth in the new grid')

% save grid
save([conf.Data.BaseDir '/grid.mat'],'Lat','Lon','BottomDepth','mask_water');
fprintf('saving grid\n');

timestamps = [];

total_outfiles = conf.OGCM.step * nFiles;
count = conf.OGCM.step * nFiles;


%% iterating through files
for n=1:nFiles

	ncfile =fnames{n};

	%Depth = double(ncread(ncfile, 'depth'));
	u_all = double(ncread(ncfile,'u'));
	v_all = double(ncread(ncfile,'v'));
	%w_all = double(ncread(ncfile,'velocity_W'));
	ssh_all = double(ncread(ncfile,'ssh'));


	fTimestamps = double(ncread(ncfile,'time'));

	for i=1:length(fTimestamps)-1
		t = fTimestamps(i);

% Time
		time = t0+t/60/60/24; % IJR 14Nov19 the variable time should be double, not duration
		time_str = datestr(time,'dd-mmm-yyyy HH:MM:SS');
		timestamps = [timestamps, time];
		fprintf('changing format for time %s\n',time_str);
	%	alltimes_str = [alltimes_str, time_str];

		E = ssh_all(:,:,i); % ssh is meters above hydrographic zero 
		E(isnan(E))=0;
		%E = E.'; %'

	% Depth at each grid point

% 		depth = nan(numlat,numlon,numlvl-num_cartesian_layers);
% 		for lv=1:length(ZZ)
% 			depth(:,:,lv) = double(-ZZ(lv)*(cartesian_depth(:,:,1)+E)+E);
% 		end
		
		u = u_all(:,:,end,i); %first layer 030420
		v = v_all(:,:,end,i);
		%w = w_all(:,:,end,i);
		
		%w(isnan(w))=0;
		u(isnan(u))=0;
		v(isnan(v))=0;
        
%         %IJR 12Nov19 I've moved this here
% 		u = flip(u,3); % velocites have deeper depths indexed first (I think)
% 		v = flip(v,3);
% 		w = flip(w,3);

        % Cut at grid size
        u=u(pos1Lon:posEndLon,pos1Lat:posEndLat,:);
        v=v(pos1Lon:posEndLon,pos1Lat:posEndLat,:);
        %w=w(pos1Lon:posEndLon,pos1Lat:posEndLat,:);
        E=E(pos1Lon:posEndLon,pos1Lat:posEndLat);

		% The velocity at the new first layer depth=0 equal to the velocity at
    % the old first sigma layer
       

		u = permute(u,[2 1]);
		v = permute(v,[2 1]);
		%w = permute(w,[2 1]);
        E = permute(E,[2 1]);


		 % From m/s to cm/s
    
   		u=u*100;
   	 	v=v*100;
    	%w=w*100;

        
%     	depth = cat(3,depth,-cartesian_depth(:,:,2:end)); %IJR14Nov19 cartasian_depth must be negative
%                                                           % IJR 030420 from the second layer
%         depth = depth-repmat(depth(:,:,1),[1,1,size(depth,3)]); %Change in the reference system (surface constant, varying bottom)
        
    	file_index = total_outfiles + 1 - count;

    	save([conf.Data.BaseDir '/TrackMPDInput' num2str(file_index) '.mat'],'u','v','E','time','time_str');
    	count = count -1;

	end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Read and save Time information    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save([conf.Data.BaseDir '/timestamps.mat'],'timestamps');
fprintf('saving timestamps\n');





