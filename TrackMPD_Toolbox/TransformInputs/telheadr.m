function m = telheadr(FILENAME)
%****M* Telemac/telheadr.m
% NAME
% telheadr.m
% 
% PURPOSE
% Read Telemac results file header (Seraphin and Leonard format).
% Works with both 2D and 3D files.
%
% USAGE
%       m = telheadr(FILENAME)
% 
% INPUTS
%       FILENAME = a filename (or file id number) for seraphin or leonard
%       file. If no inputs are given a menu will pop up asking for a file
%
% RESULTS
%       The output (m) is a structure array  containing the following
%       info (for example):
%    
%      filename: 'd:\seraphin_example.res'  -> filename
%           fid: 6                          -> the file ID
%         title: [1x80 char]                -> model run title
%           NBV: 4                          -> number of variables in file
%          RECV: {4x1 cell}                 -> variables names (cell array)
%          type: 'seraphin'                 -> file type (seraphin or leonard)
%        IPARAM: [10x1 double]              -> logical flag array
%         IDATE: []                         -> optional start date [y m d H M S]
%         NPLAN: 1                          -> number of vertical layers (2D = 1)
%         NELEM: 23071                      -> number of elements in mesh
%         NPOIN: 12304                      -> number of nodes in mesh
%           NDP: 3                          -> number of points per element
%          IKLE: [23071x3 int32]            -> interconnectivity array (NELEM * 3)
%         IPOBO: [12304x1 double]           -> boundary numbering (NPOIN * 1)
%           XYZ: [12304x2 double]           -> X and Y coorinates (NPOIN *2)
%         first: 1                          -> has a value of 1 if it is the first timestep
%          last: 0                          -> has a value of 1 if it is the last timestep
%        RESULT: [12304x4 double]           -> The timestep results array (NPOIN * NBV)
%     startfpos: 424868                     -> file position of first timestep (bytes) 
%       len1rec: 196908                     -> file length of each timestep (bytes)
%      timestep: 0                          -> the current timestep
%        NSTEPS: 50                         -> total number of timesteps in the file
%            DT: 1800                       -> the timestep interval (in seconds)
%
% EXAMPLE 
%       In this example, the 3rd and 5th timesteps are extracted from a telemac results
%       file and rewritten out to a new file.
%
%       m = telheadr('seraphin_example.res');
%       fid = telheadw(m,'extracted_steps.sel');
%       for i=[3 5]
%           m = telstepr(m,i);
%           fid = telstepw(m,fid);
%       end
%       fclose(m.fid);
%       fclose(fid);
%
% NOTES
%
% SEE ALSO
% telstepr.m telheadw.m telstepw.m
%
% AUTHOR
% Thomas Benson
% HR Wallingford
% t.benson@hrwallingford.co.uk
%
% DATE
% 13-Aug-2009
%
%****

if nargin<1
   FILENAME = [];
end


m=[];

if isempty(FILENAME)
    % ask for results filename
    [res_file, pth]=uigetfile('*.*','Select Telemac results file');
    if res_file==0, return, end % user aborted
    FILENAME = [pth res_file];
end

% if FILENAME is a string, assume it is a filename and read it directly, otherwise assume it is a
%file pointer (big endian)
if ischar(FILENAME)
    m.filename = FILENAME;
    fid = fopen(FILENAME,'r','b');
end

% add fid to the structure
m.fid = fid;

% Title
rectag1 = fread(fid,1,'int32');  %get first record start tag
m.title = (fread(fid,rectag1,'*char'))';
rectag1 = fread(fid,1,'int32');  %get first record end tag

% Number of 2D result variables (NBV)
rectag2 = fread(fid,1,'int32');  %get second record start tag
m.NBV = fread(fid,1,'int32');
NBVX = fread(fid,1,'int32');
% if NBVX ~= 0
%     error('Second mesh discretisation exists: this is not supported')
% end
rectag2 = fread(fid,1,'int32');  %get second record end tag
% fprintf('\nThere are %d variables in the TELEMAC results file, named:\n',NBV);

% Name of each variable
m.RECV = cell(m.NBV,1);
for i = 1:m.NBV
    rectag3 = fread(fid,1,'int32');  %get third record start tag
    m.RECV{i} = fread(fid,rectag3,'*char')';
    rectag3 = fread(fid,1,'int32');  %get third record end tag
end
% fprintf('%s\n',RECV{:}); % display the 2D variable names

% Parameter array
rectag4 = fread(fid,1,'int32');  %get fourth record start tag

if rectag4==40 % SERAPHIN FORMAT FILE

    disp('Seraphin format file')
    
    % N.B. the 3D format files contain elevations for nodes as the first
    % variable
    m.type = 'seraphin';

    m.IPARAM = zeros(10,1);
    m.IPARAM = fread(fid,10,'int32');
    rectag4 = fread(fid,1,'int32');  %get fourth record end tag

   % Date
    if m.IPARAM(10) == 1
        % extra date record will be here
        % 6 element data vector (Y M D H M S)
        rectag5 = fread(fid,1,'int32');  %#ok<NASGU> %get fifth record end tag
        m.IDATE = zeros(1,6);
        m.IDATE(:) = fread(fid,6,'int32');
        rectag5 = fread(fid,1,'int32');  %#ok<NASGU> %get fifth record end tag
    else
        m.IDATE = [];
    end

    % number of layers (check for TELEMAC3D output)
    if m.IPARAM(7) > 0
        m.NPLAN = m.IPARAM(7);
    else
        m.NPLAN = 1;
    end


    % Number of 2D elements and nodes
    rectag6 = fread(fid,1,'int32');  %get sixth record start tag
    m.NELEM = fread(fid,1,'int32');
    m.NPOIN = fread(fid,1,'int32');
    m.NDP = fread(fid,1,'int32');
    not_used = fread(fid,1,'int32');
    rectag6 = fread(fid,1,'int32');  %get sixth record end tag;

    % fprintf('\nNumber of elements = %d\n',NELEM);
    % fprintf('Number of nodes = %d\n',NPOIN);

    % pre-dimension some more arrays
    %if m.NPLAN>1
    %   m.IKLE = repmat(int32(0),m.NELEM,m.NDP);
    %else
    m.IKLE = repmat(int32(0),m.NDP,m.NELEM);
    %end

    m.IPOBO = zeros(m.NPOIN,1);
    m.XYZ = zeros(m.NPOIN,2);

    % Connectivity table
    rectag7 = fread(fid,1,'int32');  %get sixth record start tag
    % now read in the IKLE array, which will be of dimension (NDP,NELEM)
    % check that tag/4 = NDP*NELEM
    if rectag7/4 ~= m.NDP*m.NELEM
        error('Error: inconsistent dimension of IKLE array')
    else
        %    if m.NPLAN==1  % 2D file
        % read the array
        m.IKLE(:) = fread(fid,[m.NDP m.NELEM],'int32');
        m.IKLE=m.IKLE';   % transpose it
        %    else % 3D file
        %        m.IKLE(:) = fread(fid,[m.NELEM m.NDP],'int32');
        %    end

    end
    rectag7 = fread(fid,1,'int32');  %get sixth record end tag

    % Boundary nodes
    rectag8 = fread(fid,1,'int32');  %get 7th record start tag
    % check that tag/4 = NPOIN
    if rectag8/4 ~= m.NPOIN
        error('Error: inconsistent dimension of IPOBO array')
    else
        % read the array
        m.IPOBO(:) = fread(fid,m.NPOIN,'int32');
    end
    rectag8 = fread(fid,1,'int32');  %get 7th record end tag

    % Node X coords
    rectag9 = fread(fid,1,'int32');  %get 8th record start tag
    % check that tag/4 = NPOIN
    if rectag9/4 ~= m.NPOIN
        error('Error: inconsistent dimension of 2D X abscissa array')
    else
        % read the array, where i = 1:NPOIN
        m.XYZ(:,1) = fread(fid,m.NPOIN,'float32');
    end
    rectag9 = fread(fid,1,'int32');  %get 8th record end tag

    % Node Y coords
    rectag10 = fread(fid,1,'int32');  %get 9th record start tag
    % check that tag/4 = NPOIN
    if rectag10/4 ~= m.NPOIN
        error('Error: inconsistent dimension of 2D Y ordinate array')
    else
        % read the array, where i = 1:NPOIN
        m.XYZ(:,2) = fread(fid,m.NPOIN,'float32');
    end
    rectag10 = fread(fid,1,'int32');  %get 9th record end tag


    if m.NPLAN==1
        fprintf('\nThere are %d variables in the 2D results file, named:\n',m.NBV);
    else
        fprintf('\nThere are %d variables in the 3D results file, named:\n',m.NBV);
    end
    fprintf('%s\n',m.RECV{:}); % display the variable names
    fprintf('\nNumber of elements = %d\n',m.NELEM);
    fprintf('Number of nodes = %d\n',m.NPOIN);
    fprintf('Number of vertical planes = %d\n',m.NPLAN);
    
elseif rectag4==20 % LEONARD FORMAT
    
    % N.B. the 3D format files contain elevations for nodes as the first
    % variable
    m.type = 'leonard';
    
    m.MESHSIZE = zeros(5,1);
    m.MESHSIZE = fread(fid,5,'int32');
    m.NPOIN = m.MESHSIZE(1)*m.MESHSIZE(2);
    rectag4 = fread(fid,1,'int32');  %get fourth record end tag

    % tokens (only 2nd one used)
    rectag5 = fread(fid,1,'int32');  %get fourth record start tag
    m.IPARAM = zeros(10,1);
    m.IPARAM = fread(fid,10,'int32');
    rectag5 = fread(fid,1,'int32');  %get fourth record end tag

    % X coordinates
    rectag6 = fread(fid,1,'int32');  %get fifth record start tag
    m.X = zeros(m.MESHSIZE(1),m.MESHSIZE(2));
    m.X(:) = fread(fid,m.NPOIN,'float32');
    m.X = m.X';
    rectag6 = fread(fid,1,'int32');  %get fifth record end tag;

    % Y coordinates
    rectag7 = fread(fid,1,'int32');  %get fifth record start tag
    m.Y = zeros(m.MESHSIZE(1),m.MESHSIZE(2));
    m.Y(:) = fread(fid,m.NPOIN,'float32');
    m.Y = m.Y';
    rectag7 = fread(fid,1,'int32');  %get fifth record end tag;

    % Inside domain mask
    rectag8 = fread(fid,1,'int32');  %get fifth record start tag
    m.INDIC = zeros(m.MESHSIZE(1),m.MESHSIZE(2));
    m.INDIC(:) = fread(fid,m.NPOIN,'int32');
    m.INDIC = m.INDIC';
    rectag8 = fread(fid,1,'int32');  %get fifth record end tag;


    fprintf('\nThere are %d variables in the results file, named:\n',m.NBV);
    fprintf('%s\n',m.RECV{:}); % display the 2D variable names
    %fprintf('\nNumber of elements = %d\n',m.NELEM);
    fprintf('Number of nodes = %d\n',m.NPOIN);


else
    error('Unknown telemac file format')
end

% add two tokens for the first and last variable
m.first = 1;
m.last  = 0;

% % also try to work out number of timesteps from the length of the file
% read first timstep to get length of 1 record
m = telstepr(m);

% forward wind to end of file
fseek(m.fid,0,'eof');
fpos = ftell(m.fid);

% calc the size
try m.len1rec;
    m.NSTEPS = floor((fpos-m.startfpos)/m.len1rec);
catch
    m.NSTEPS=0;
end
    
fprintf('Number of timesteps = %d\n',m.NSTEPS);

% rewind to start of first timestep record
%fseek(m.fid,m.startfpos,'bof');

% if we have at least 2 timesteps, calculate the timestep interval
if m.NSTEPS>1
    
    m = telstepr(m,1);
    T0 = m.AT;
    m = telstepr(m,2);
    T1 = m.AT;

    % calculate the timestep
    m.DT = T1-T0;

else
    
    m.DT = 0;
    
end

% rewind to start of first timestep record
fseek(m.fid,m.startfpos,'bof');
    
% restate the two tokens for the first and last variable
m.first = 1;
m.last  = 0;
m.timestep = 0;
try m.AT;
    m=rmfield(m,'AT');
end
