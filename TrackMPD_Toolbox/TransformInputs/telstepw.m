function fid = telstepw(m,fid)
%****M* Telemac/telstepw.m
% NAME
% telstepw.m
% 
% PURPOSE
% Write Telemac results file header (Seraphin or Leonard format).
% Works with both 2D and 3D files.
%
% USAGE
%       fid = telstepw(m,fid)
% 
% INPUTS
%       m = a structure array as produced by telheadr.m (see telheadr.m)
%
%       fid = usually a file pointer to a serpahin or leonard
%       file (or can be a FILENAME).
%
% RESULTS
%       fid = the updated file pointer to the results file
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
% (1) You will need to run telheadw.m before this.
% (2) It is important that the time variable (m.AT) is updated each time before
% running this program. i.e. each timestep must have a unique time which is
% increasing! 
% (3) If no output is defined, the file will be closed. 
%
% SEE ALSO
% telstepr.m telheadr.m telheadw.m
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

% write telemac timestep info to seraphin or leonard file
% note: run telemacleaderwrite.m before running this function.
%
% m = telemac struct array (seraphin or leonard format)
%            containing timestep info
% fid = seraphin filename to write data to, OR file id number
%
% Written by:
% Tom Benson 
% HRwallingford
% Nov-2006

if nargin<2

    % here, the fid within the structure array is used
%     try m.fid; % check it's there
%     catch % if not, ask for a filename (MUST CONTAIN LEADER INFO ALREADY!)
        % ask for results filename
        [res_file, pth]=uigetfile('*.*','Select Telemac leader file for appending step data');
        if res_file==0, return, end % user aborted
        fid = fopen([pth res_file],'ab','b');
%     end

else

    % if a file or fid has been input manually
    if ischar(fid) % NOTE: if fid is a name, it must already contain leader data
        fid = fopen(fid,'ab','b');
    %else
        %fid = FILENAME;
    end

end
    
if strcmp(m.type,'seraphin')
    
    % write time for this timestep
    rectag10 = 4;
    fwrite(fid,rectag10,'int32');  %set 10th record start tag
    fwrite(fid,m.AT,'float32');
    fwrite(fid,rectag10,'int32');  %set 10th record end tag

    % write RESULT for this timestep
    rectag11 = m.NPOIN * 4;
    for i=1:m.NBV
        fwrite(fid,rectag11,'int32');  %set 11th record start tag
        fwrite(fid,m.RESULT(:,i),'float32');
        fwrite(fid,rectag11,'int32');  %set 11th record end tag
    end

elseif strcmp(m.type,'leonard')
   
    % write time for this timestep
    rectag10 = 4;
    fwrite(fid,rectag10,'int32');  %set 10th record start tag
    fwrite(fid,m.AT,'float32');
    fwrite(fid,rectag10,'int32');  %set 10th record end tag
    
    rectag11 = m.MESHSIZE(1) * m.MESHSIZE(2) * 4;
    % write RESULT for this timestep
    for i = 1:m.NBV
        fwrite(fid,rectag11,'int32');  %set 11th record start tag
        fwrite(fid,m.RESULT(:,:,i),'float32');
        fwrite(fid,rectag11,'int32');  %set 11th record start tag
    end
    
else
    error('Unrecognised Telemac file format')
    
end

if nargout==0
   fclose(fid);
end