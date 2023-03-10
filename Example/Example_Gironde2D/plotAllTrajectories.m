clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%% Modify %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
path_data='Outputs';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*


FileList = dir([path_data '/Gironde*.mat']);


for i=1:length(FileList)
  filename = FileList(i).name;
  
  %% Open data (Trajectories)
  load([path_data '/' filename]);
  
  Lon=TRAJ.Lon;
  Lat=TRAJ.Lat;
  Depth=TRAJ.Depth;
  FinalLonLatDepth=TRAJ.FinalLonLatDepth;
  Dt = (TRAJ.TimeStamp(3)-TRAJ.TimeStamp(2))*3600*24;
  U = (Lon(:,2:end)-Lon(:,1:end-1))./Dt;
  V = (Lat(:,2:end)-Lat(:,1:end-1))./Dt;
  Vel = zeros(size(Lon));
  Vel(:,1:end-1) = sqrt(U.^2+V.^2);
  Vel(:,end)=Vel(:,end-1);
  
  
  % Number of particles and of time stamps
  [NumParticles TimeSteps]=size(Lat);
  
  
  %% Plot trajectory
  
  % figure;
  figure('units','normalized','outerposition',[0.1+0.05*(i-1) 0.1 0.5 0.6]);
  
  %plot domain
  fill(TRAJ.Domain.x,TRAJ.Domain.y,[0.8 0.8 0.8]);
  hold on;
  Opt.Marker = 'o';
  Opt.MarkerSize = 6;
  Opt.ColorMap = jet(128);
  Opt.CRange = [0 2];
  Opt.Line = 1;
  Opt.LineWidth = 2;
  %plot trajectories
  for i=1:NumParticles
%    plot(Lon(i,:),Lat(i,:),'-x','LineWidth',1); hold on; %Trajectory
    scatter(Lon(i,:),Lat(i,:),[],Vel(i,:),'filled');
    colorbar
    hold on
    plot(Lon(i,1),Lat(i,1),'x','MarkerFaceColor','b','LineWidth',2,'MarkerSize',7,'MarkerEdgeColor','k'); %Release point
    plot(FinalLonLatDepth(i,1),FinalLonLatDepth(i,2),'.','color','r','MarkerSize',10); %Fate
  end
  %title(filename)
  set(gca,'FontSize',16)
  
  % ylim([6494000 6514000]); xlim([382000 388000]); %Embouchure
  % ylim([6450000 6550000]); xlim([350000 420000]); %Estuary
  ylim([6430000 6520000]); xlim([360000 440000]);
  axis equal
  
  title(strrep(filename,"_"," "))
  
end