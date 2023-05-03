%Script to check that the domain file is correct
% I.Jalón Rojas April 2023

domain=load('domain.txt');

LON=domain(:,1)';
LAT=domain(:,2)';

aux=find(isnan(LON)==1);
aux=[0 aux];

figure;
for i=1:length(aux)-1
    fill(LON(aux(i)+1:aux(i+1)-1),LAT(aux(i)+1:aux(i+1)-1),[0.5 0.5 0.5])
    hold on;
    
end

%check if land is plotted by closed grey polygons (check TrackMPD tutorial
%on GitHub for more details)
