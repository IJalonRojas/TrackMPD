function fileList = getAllFiles(dirName, fileExtension, appendFullPath)
% GETALLFILES  lists files and folders in the specific folder
%
% Usage: fileList = getAllFiles(dirName, fileExtension, appendFullPath)
% 
% Inputs
% ------
% dirName:      Name of the directory
% fileExtension: File extension
% appendFullPath: append the Full Path
%
% Outputs
% -------
% fileList: Files name
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	$Id: getAllFiles.m 001 2018-03-05 18:53:15Z ef $
%
% Copyright (C) 2014-2018 Erick Fredj
% Licence: GPL (Gnu Public License)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  dirData = dir([dirName filesep fileExtension]);      %# Get the data for the current directory
  dirWithSubFolders = dir(dirName);
  dirIndex = [dirWithSubFolders.isdir];  %# Find the index for directories
  fileList = {dirData.name}';  %'# Get a list of the files
  if ~isempty(fileList)
    if appendFullPath
      fileList = cellfun(@(x) fullfile(dirName,x),...  %# Prepend path to files
                       fileList,'UniformOutput',false);
    end
  end
  subDirs = {dirWithSubFolders(dirIndex).name};  %# Get a list of the subdirectories
  validIndex = ~ismember(subDirs,{'.','..'});  %# Find index of subdirectories
                                               %#   that are not '.' or '..'
  for iDir = find(validIndex)                  %# Loop over valid subdirectories
    nextDir = fullfile(dirName,subDirs{iDir});    %# Get the subdirectory path
    fileList = [fileList; getAllFiles(nextDir, fileExtension, appendFullPath)];  %# Recursively call getAllFiles
  end

end