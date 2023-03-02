function h=telplot(INFILE,VAR,TIMESTEP,LIMS)

% Plots a patch for a TELEMAC2D triangular, finite element style mesh
% (i.e. a Seraphin format file)
%
% usage:
%  h=telplot(INFILE,VAR,TIMESTEP,LIMS)
%
% where:
% INFILE = filename of Seraphin file (setting to [], or no inputs, will invoke a uigetfile menu)
% VAR =  variable number (default = 1)
% TIMESTEP = timestep number (default = 1)
% LIMS  = [min max] limits for colour bar, (default are limits of data)
%
% The function returns the handle to the patch object
%
% example:
%  telplot('seraphin_example.res',4,1)
%
% author: Thomas Benson, HR Wallingford
% email: t.benson@hrwallingford.co.uk
% release date: 13-Aug-2009

% check inputs
if nargin<1
    INFILE = [];
end
if nargin<2
    VAR=[];
end
if nargin<3
    TIMESTEP=[];
end
if nargin<4
    LIMS=[];
end
if isempty(VAR)
    VAR=1;
end
if isempty(TIMESTEP)
    TIMESTEP=1;
end

% load the header
m=telheadr(INFILE);

% load the timestep
m=telstepr(m,TIMESTEP);

% check inputs are in range
if VAR>m.NBV
    error(['VAR (input 2) is too large. There are only ' num2str(m.NBV) ' variables in the file.'])
end
if TIMESTEP>m.NSTEPS
    error(['TIMESTEP (input 3) is too large. There are only ' num2str(m.NSTEPS) ' timesteps in the file.'])
end

% plot the patch
h=patch('faces',m.IKLE,'vertices',m.XYZ,'FaceVertexCData',m.RESULT(:,VAR), ...
    'FaceColor','interp','EdgeColor','none','linewidth',0.01);
axis equal
axis tight
set(gca,'box','on')
xlabel('East (m)')
ylabel('North (m)')
colorbar

% change the colour scaling
if isempty(LIMS)
    caxis([min(m.RESULT(:,VAR)) max(m.RESULT(:,VAR))]);
else
    caxis(LIMS);
end

% close the file
fclose(m.fid);
