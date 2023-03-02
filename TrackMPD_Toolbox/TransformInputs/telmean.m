function m=telmean(INFILE,OUTFILE)
% Example function that reads a Telemac2D (Seraphin) result file
% and averages the data for all the timesteps, then outputs to a new 
% Seraphin file
%
% usage: telmean(INFILE,OUTFILE)
%
% INFILE = input filename (optional - a uimenu will appear if no inputs)
% OUTFILE = output filename (optional - an output file will be created if omitted)
%
% other functions that are called:
%
% telheadr.m - for reading Seraphin file header
% telheadw.m - for reading Seraphin file timestep
% telstepr.m - for writing Seraphin file header
% telstepw.m - for writing Seraphin file timestep
%
% author: Thomas Benson, HR Wallingford
% email: t.benson@hrwallingford.co.uk
% release date: 13-Aug-2009

% check inputs inputs
if nargin<1
    INFILE=[];
end
if nargin<2
    OUTFILE = [];
end

% read the file header
m = telheadr(INFILE);

% get the filename parts
INFILE = m.filename;
[PATHSTR,NAME,EXT] = fileparts(INFILE);

% create an averaging array
AVGARRAY = repmat(0,size((m.RESULT)));

fprintf('\nNumber of timesteps: %d\n',m.NSTEPS)
fprintf('Current timestep: %10d',0)

% loop through the timesteps
for i=1:m.NSTEPS
    
    fprintf('\b\b\b\b\b\b\b\b\b\b%10d',i);
    
    % read the timestep
    m = telstepr(m,i);

    % add the data to the avergaing array
    AVGARRAY = m.RESULT+AVGARRAY;

end
fprintf('\n');

% make sure the input file is closed
fclose(m.fid);

% divide by the number of timesteps to get the averages
% and put back into the structure result array
m.RESULT = AVGARRAY./m.NSTEPS;

% reset the time to zero for outputting
m.AT = 0;

% if no output file is specified, create a filename
if isempty(OUTFILE)
    OUTFILE = [NAME '_mean' EXT];
end

% if the function is called without an output, write to file.
if nargout<1
    % write the header to output
    fid = telheadw(m,OUTFILE);
    % write the timestep (averaged data) to output
    fid = telstepw(m,fid);
    % close the file
    fclose(fid);
end

