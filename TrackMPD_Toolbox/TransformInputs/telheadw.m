function fid = telheadw(m,FILENAME)
%****M* Telemac/telheadw.m
% NAME
% telheadw.m
% 
% PURPOSE
% Write Telemac results file header (Seraphin and Leonard format).
% Works with both 2D and 3D files.
%
% USAGE
%       fid = telheadw(m,FILENAME)
% 
% INPUTS
%       m = a structure array as produced by telheadr.m
%
%       FILENAME = a filename (or file id number) for serpahin or leonard
%       file.
%
% RESULTS
%       fid = a file pointer to the written header
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
% telstepr.m telheadr.m telstepw.m
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
    
    disp('Writing seraphin format leader information')
    
    % write out a seraphin file with defined variables
    % the m array is a structure containing (for example):
    %            title: [1x80 char]
    %             NBV: 5
    %            RECV: {5x1 cell}
    %           IPARAM: [10x1 double]
    %            IDATE: []
    %            NPLAN: 1
    %           NELEM: 22265
    %           NPOIN: 11562
    %              NDP: 3
    %            IKLE: [22265x3 int32]
    %            IPOBO: [11562x1 double]
    %             XYZ: [11562x2 double]
    %

    % check if this is a 3D file (IKLE will be NPOIN * 6)
    try
        m.IKLE;  % see if we are using the standard convention
    catch
        try
            m.IKLE3; % contains 3D type names
            % swap to normal convention
            m.IKLE = m.IKLE3;
            m.NBV  = m.NBV3;
            m.RECV = m.RECV3;
            m.IPARAM = m.IPARAM3;
            m.NPOIN = m.NPOIN3;
            m.NELEM = m.NELEM3;
            %m.fid = m.fid3;

        catch
            error('Structure array does not contain required info')
        end
    end

    if size(m.IKLE,2)==6
        disp('3D file format detected')
        try
            m.NPLAN;
            if m.IPARAM(7)~=m.NPLAN
                error('NPLAN and IPARAM(7) are inconsistent...check number of levels')
            end
        catch
            if m.IPARAM(7)<=1
                error('Need to know the number of levels in IPARAM(7)')
            end
            %m.NPLAN = m.IPARAM(7);
        end
    end

    fid = fopen(FILENAME,'wb','b');
    %fid = fopen(FILENAME,'r+','b');
    rectag1 = length(m.title);
    if rectag1 > 80	% restrict record length to 80 chars
        rectag1 = 80;
        %message = 'Study title truncated to: '
        m.title = sscanf(m.title,'%c',80);
    elseif rectag1 < 80 % pack out to 80 characters
        paddingreqd = 80-rectag1;
        pad = 32 .* ones(1,paddingreqd);
        pad = char(pad); % make up some spaces
        m.title = [m.title pad];
        rectag1 = 80;
    end

    % write title string to first record
    fwrite(fid,rectag1,'int32');  %set 1st record start tag
    fwrite(fid,m.title,'char')'; %write title string
    fwrite(fid,rectag1,'int32');  %set 1st record end tag

    % NBV records (add to the existing variables)
    %m.NBV = size(VAR,2);

    rectag2 = 8; % since they are written as two 4 byte integers
    fwrite(fid,rectag2,'int32');  %set 2nd record start tag
    fwrite(fid,m.NBV,'int32');
    fwrite(fid,0,'int32');
    fwrite(fid,rectag2,'int32');  %set 2nd record end tag

    % now write variable name strings for NBV
    for i=1:length(m.RECV)
        rectag3 = length(m.RECV{i});
        if rectag3 > 32	% restrict record length to 32 chars
            rectag3 = 32;
            %message = 'Variable name truncated to: '
            m.RECV{i} = sscanf(m.RECV{i},'%c',32);
        else % pack to 32 characters
            paddingreqd = 32-rectag3;
            pad = 32 .* ones(1,paddingreqd);
            pad = char(pad); % make up some spaces
            %message = 'Variable name padded to: '
            m.RECV{i} = [m.RECV{i} pad];
            rectag3 = 32;
        end
    end

    %m.RECV = [m.RECV;EXTRA_RECV2];
    %rectag3 = 32;
    for i = 1:length(m.RECV)
        % write out variable name and units to 3rd record
        fwrite(fid,rectag3,'int32');  %set 3rd record start tag
        fwrite(fid,m.RECV{i},'char')';
        fwrite(fid,rectag3,'int32');  %set 3rd record end tag
    end

    % now write out an integer vector, IPARAM
    % there is no equivalent for this in SMS/RMA
    % but the values are prescribed anyway
    rectag4 = 40;   % ten 4 byte integers
    %IPARAM = [1 0 0 0 0 0 0 0 0 0];
    fwrite(fid,rectag4,'int32');  %set 4th record start tag
    count = fwrite(fid,m.IPARAM,'int32')';
    fwrite(fid,rectag4,'int32');  %set 4th record end tag

    % now we need to find out number of elements, points and points per element
    % NELEM = number of elements
    % NPOIN = number of points
    % NDP = number of points per element (3 for triangles)
    % unused = 1 (an unused integer value)
    
    % Date
    if isfield(m,'IDATE') && ~isempty(m.IDATE)%m.IPARAM(10) == 1
        m.IPARAM(10)=1; % make sure this is set to 1
        % extra date record
        % 6 element data vector (Y M D H M S)
        rectagextra = 24;
        fwrite(fid,rectagextra,'int32'); 
        fwrite(fid,m.IDATE,'int32');
        fwrite(fid,rectagextra,'int32'); 
    else
        m.IDATE = [];
    end
    
    unused = 1;

    % now write out the summary field
    rectag5 = 16; % four 4 byte integers
    fwrite(fid,rectag5,'int32');  %set 5th record start tag
    fwrite(fid,m.NELEM,'int32');
    fwrite(fid,m.NPOIN,'int32');
    fwrite(fid,m.NDP,'int32');
    fwrite(fid,unused,'int32');
    fwrite(fid,rectag5,'int32');  %set 5th record end tag


    % now write the IKLE table (element connectivity)
    % first compute record length, which will be NDP * NELEM * 4 bytes
    rectag6 = m.NDP * m.NELEM * 4;
    fwrite(fid,rectag6,'int32');  %set 6th record start tag
    j = 1:m.NDP;
    i = 1:m.NELEM;
    m.IKLE = m.IKLE';
    fwrite(fid,m.IKLE(j,i),'int32');
    fwrite(fid,rectag6,'int32');  %set 6th record end tag

    % record length is NPOIN * 4
    rectag7 = m.NPOIN * 4;
    fwrite(fid,rectag7,'int32');  %set 7th record start tag
    fwrite(fid,m.IPOBO,'int32');
    fwrite(fid,rectag7,'int32');  %set 7th record end tag

    % write X and Y values
    rectag8 = m.NPOIN * 4;
    fwrite(fid,rectag8,'int32');  %set 8th record start tag
    fwrite(fid,m.XYZ(:,1),'float32');
    fwrite(fid,rectag8,'int32');  %set 8th record end tag
    rectag9 = rectag8;
    fwrite(fid,rectag9,'int32');  %set 8th record start tag
    fwrite(fid,m.XYZ(:,2),'float32');
    fwrite(fid,rectag9,'int32');  %set 8th record end tag


elseif strcmp(m.type,'leonard')

       
    disp('Writing leonard format leader information')
    
    fid = fopen(FILENAME,'wb','b');

    rectag1 = length(m.title);
    if rectag1 > 80	% restrict record length to 80 chars
        rectag1 = 80;
        %message = 'Study title truncated to: '
        m.title = sscanf(m.title,'%c',80);
    elseif rectag1 < 80 % pack out to 80 characters
        paddingreqd = 80-rectag1;
        pad = 32 .* ones(1,paddingreqd);
        pad = char(pad); % make up some spaces
        m.title = [m.title pad];
        rectag1 = 80;
    end
    % write title string to first record
    fwrite(fid,rectag1,'int32');  %set 1st record start tag
    fwrite(fid,m.title,'char')'; %write title string
    fwrite(fid,rectag1,'int32');  %set 1st record end tag

    % NBV records (add to the existing variables)
    %m.NBV = size(VAR,2);

    rectag2 = 8; % since they are written as two 4 byte integers
    fwrite(fid,rectag2,'int32');  %set 2nd record start tag
    fwrite(fid,m.NBV,'int32');
    fwrite(fid,0,'int32');
    fwrite(fid,rectag2,'int32');  %set 2nd record end tag

    % now write variable name strings for NBV
    for i=1:length(m.RECV)
        rectag3 = length(m.RECV{i});
        if rectag3 > 32	% restrict record length to 32 chars
            rectag3 = 32;
            %message = 'Variable name truncated to: '
            m.RECV{i} = sscanf(m.RECV{i},'%c',32);
        else % pack to 32 characters
            paddingreqd = 32-rectag3;
            pad = 32 .* ones(1,paddingreqd);
            pad = char(pad); % make up some spaces
            %message = 'Variable name padded to: '
            m.RECV{i} = [m.RECV{i} pad];
            rectag3 = 32;
        end
    end

    %m.RECV = [m.RECV;EXTRA_RECV2];
    %rectag3 = 32;
    for i = 1:length(m.RECV)
        % write out variable name and units to 3rd record
        fwrite(fid,rectag3,'int32');  %set 3rd record start tag
        fwrite(fid,m.RECV{i},'char')';
        fwrite(fid,rectag3,'int32');  %set 3rd record end tag
    end


    % number of columns and lines in mesh
    rectag4 = 20;  %fourth record start tag 5*4 byte integers
    fwrite(fid,rectag4,'int32');
    fwrite(fid,m.MESHSIZE,'int32');
    fwrite(fid,rectag4,'int32');  %get fourth record end tag

    % tokens (only 2nd one used)
    rectag5 = 40;  %5th record start tag 10*4 byte integers
    fwrite(fid,rectag5,'int32');
    fwrite(fid,m.IPARAM,'int32');
    fwrite(fid,rectag4,'int32');

    % X coordinates
    rectag6 = m.MESHSIZE(1)*m.MESHSIZE(2)*4;
    fwrite(fid,rectag6,'int32');
    fwrite(fid,m.X','float32');
    fwrite(fid,rectag6,'int32');

    % X coordinates
    rectag7 = m.MESHSIZE(1)*m.MESHSIZE(2)*4;
    fwrite(fid,rectag7,'int32');
    fwrite(fid,m.Y','float32');
    fwrite(fid,rectag7,'int32');

    % Inside domain mask
    rectag8 = m.MESHSIZE(1)*m.MESHSIZE(2)*4;
    fwrite(fid,rectag8,'int32');
    fwrite(fid,m.INDIC','int32')
    fwrite(fid,rectag8,'int32');

end

if nargout==0
    fclose(fid);
end