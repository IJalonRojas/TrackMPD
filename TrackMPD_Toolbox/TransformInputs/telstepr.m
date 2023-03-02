function m = telstepr(m,TIMESTEP)
%****M* Telemac/telstepr.m
% NAME
% telstepr.m
% 
% PURPOSE
% Read Telemac results file timestep information (Seraphin and Leonard format).
% Works with both 2D and 3D files.
%
% USAGE
%       m = telstepr(m,TIMESTEP)
% 
% INPUTS
%       m = a structure array as produced by telheadr.m
%       TIMESTEP = the number of the timestep to read
%
% RESULTS
%       The program will update the RESULTS and time (AT) arrays within the
%       structure array m with the information held in the Telemac file at
%       the specified timestep number (TIMESTEP)
%
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
% (1) the first call to this function may add variables to the 
% m array which will be used on subsequent calls for speeding
% things up. (i.e. length of 1 record, file position of first record , 
% last timestep read,etc.)
% (2) Once you have made one call to the file using telheadr.m you can
% make multiple calls using telstepr.m without re-running telheadr. To
% close the file, use 'fclose all', or 'fclose(m.fid)'.
%
% SEE ALSO
% telheadr.m telheadw.m telstepw.m
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

% check for file type, and assign automatically by looking for INDIC
try m.type;
catch
    try m.INDIC;
        m.type = 'leonard';
    catch
        try m.NELEM;
            m.type = 'seraphin';
        catch
            error('Structure array does not contain recognised data format')
        end
    end
end

if strcmp(m.type,'seraphin')

   
    if m.first % if this is the first time we are reading the step info

        % prepare results array
        m.RESULT = zeros(m.NPOIN,m.NBV);

        % get the file location
        m.startfpos = ftell(m.fid);
    end

    if m.first||nargin<2 % if this is first one, or if we just want to read next timestep

        % get first 2D time start tag
        steptag = fread(m.fid,1,'int32');

        if isempty(steptag)
            if m.first
                warning('No results data in file')
            end
            m.last = 1;
            return
        end

        % Time info
        m.AT = fread(m.fid,1,'float32');     % read time (in seconds)

        % get 2D time end tag
        tag = fread(m.fid,1,'int32');

        % 2D variable data
        for i = 1:m.NBV
            tag = fread(m.fid,1,'int32'); % get 2D results start tag
            if isempty(tag),
                continue,
            end;  % this may happen if reading incomplete results file
            try
                m.RESULT(:,i) = fread(m.fid,m.NPOIN,'float32'); % read data
            catch
                warning('Incomplete timestep data - returning')
                break
            end
            tag = fread(m.fid,1,'int32');  %get 2D results end tag
        end
    end

    if m.first
        % calculate record length
        m.len1rec = ftell(m.fid) - m.startfpos;
        m.first = 0; % set first to zero
    end

    % work out which timestep was read last
    m.timestep = (ftell(m.fid) - m.startfpos)/m.len1rec;
    
    % wind (fwd or rwd) to the correct position if we want a
    % specified step
    if nargin == 2
        
       
        if TIMESTEP~=m.timestep % don't need to read again if the current step

            fseek(m.fid,m.len1rec*(TIMESTEP-m.timestep-1),'cof');

            % get first 2D time start tag
            steptag = fread(m.fid,1,'int32');

            if isempty(steptag)
                m.last = 1;
                return
            end

            % find the 2D file position
            %fpos = ftell(m.fid);

            % Time
            m.AT = fread(m.fid,1,'float32');     % read time (in seconds)

            % get 2D time end tag
            tag = fread(m.fid,1,'int32');

            % 2D variable data
            for i = 1:m.NBV
                tag = fread(m.fid,1,'int32'); % get 2D results start tag
                if isempty(tag),
                    continue,
                end;  % this may happen if reading incomplete results file
                try
                    m.RESULT(:,i) = fread(m.fid,m.NPOIN,'float32'); % read data
                catch
                    warning('Incomplete timestep data - returning')
                    break
                end
                tag = fread(m.fid,1,'int32');  %get 2D results end tag
            end
        end
        % work out which timestep we have just read
        m.timestep = ((ftell(m.fid) - m.startfpos)/m.len1rec);
        
        % check for end of record (m.timestep will be less than
        % TIMESTEP)
        if m.timestep<TIMESTEP
             %error('Problem reading the correct timestep')
             fprintf('\nEnd of file reached\n')
             m.last = 1;
        end
    end


elseif strcmp(m.type,'leonard')

    if m.first % if this is the first time we are reading the step info

        % prepare results array
        m.RESULT = zeros(m.MESHSIZE(1),m.MESHSIZE(2),m.NBV);

        % get the file location
        m.startfpos = ftell(m.fid);
    end

    if m.first||nargin<2 % if this is first one, or if we just want to read next timestep

        % get first 2D time start tag
        steptag = fread(m.fid,1,'int32');

        if isempty(steptag)
            m.last = 1;
            return
        end

        % Time info
        m.AT = fread(m.fid,1,'float32');     % read time (in seconds)

        % get 2D time end tag
        tag = fread(m.fid,1,'int32');

        % 2D variable data
        for i = 1:m.NBV
            tag = fread(m.fid,1,'int32'); % get 2D results start tag
            if isempty(tag),
                continue,
            end;  % this may happen if reading incomplete results file

            m.RESULT(:,:,i) = fread(m.fid,[m.MESHSIZE(1) m.MESHSIZE(2)],'float32'); % read data

            tag = fread(m.fid,1,'int32');  %get 2D results end tag
        end

    end

    if m.first
        % calculate record length
        m.len1rec = ftell(m.fid) - m.startfpos;
        m.first = 0; % set first to zero
    end

    % work out which timestep we have just read
    m.timestep = (ftell(m.fid) - m.startfpos)/m.len1rec;

    % wind (fwd or rwd) to the correct position if we want other a
    % specified step
    if nargin == 2

        if TIMESTEP~=m.timestep % don't need to read again if the current step

            fseek(m.fid,m.len1rec*(TIMESTEP-m.timestep-1),'cof');

            % get first 2D time start tag
            steptag = fread(m.fid,1,'int32');

            if isempty(steptag)
                m.last = 1;
                return
            end

            % find the 2D file position
            %fpos = ftell(m.fid);

            % Time
            m.AT = fread(m.fid,1,'float32');     % read time (in seconds)

            % get 2D time end tag
            tag = fread(m.fid,1,'int32');

            % 2D variable data
            for i = 1:m.NBV
                tag = fread(m.fid,1,'int32'); % get 2D results start tag
                if isempty(tag),
                    continue,
                end;  % this may happen if reading incomplete results file

                m.RESULT(:,:,i) = fread(m.fid,[m.MESHSIZE(1) m.MESHSIZE(2)],'float32'); % read data

                tag = fread(m.fid,1,'int32');  %get 2D results end tag
            end
        end

        % work out which timestep we have just read (to check it is same!)
        m.timestep = (ftell(m.fid) - m.startfpos)/m.len1rec;

        % check for end of record (m.timestep will be less than
        % TIMESTEP)
        if m.timestep<TIMESTEP
             %error('Problem reading the correct timestep')
             fprintf('\nEnd of file reached\n')
             m.last = 1;
        end

    end
end















