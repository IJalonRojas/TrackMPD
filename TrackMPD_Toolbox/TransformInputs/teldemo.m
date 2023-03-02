% This is just a simple script for running the example
% files: telmean.m and telplot.m

% calculate the average of all the variables in the example result file
telmean('mersey.res','timeaverage.res')

% plot the average of variable number 3 (which is WATER DEPTH)
h=telplot('timeaverage.res',3);

% note: telplot sets the edge colour of the patch to 'none'
% so to view the mesh, set the edge colour of the elements to black 
set(h,'edgecolor','k')