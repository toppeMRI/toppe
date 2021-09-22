function [sa_full_filename] = getSAfn(sapath, varargin)
% get the full filename (including path) of the ScanArchive file located at
% sapath. The defaults assumes one ScanArchive file per folder, but will it 
% will return an individual file using the variable input 
% argument "filename_index"
%
% Inputs:
%   sapath  - string containing full path of the ScanArchive file
%
% Outputs:
%   sa_full_filename - string contraining the full filename (including
%                      path) of the ScanArchive file
%
% Optional Inputs:
%   filename_index - integer stating which file to open
%
%
% Melissa Haskell, 2021, mhask@umich.edu

%% Parse variable inputs
arg.filename_index = [];
arg = vararg_pair(arg, varargin);


%% Get files at the given path and find ScanArchive file

sapath_files = dir(sapath);  % file info at this path
name = {sapath_files.name};  % filenames
name = name(~strncmp(name, '.', 1)); % Removed files starting with '.'
name = name(contains(name, '.h5'));  % Only keep .h5 files

if numel(name) == 1
    sa_full_filename = [sapath, name{1}];
else
    if numel(arg.filename_index) == 1
        sa_full_filename = [sapath, name{arg.filename_index}];
    else
        error('Error: Folder contains multiple files, and file index not specified.')
    end
end

end