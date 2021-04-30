function [] = update_modules(tipdown_index,readout_index)
%
% Inputs:
%   tipdown_index - index of module that should be the tipdown.mod file
%   readout_index - index of module that should be the tipdown.mod file
%
%
% This script is a helper function for converting pulseq files to toppe 
% files. Specifically it:
%   1. renames the specified modules to "tipdown.mod" and "readout.mod"
%   2. makes a copy of the specified module files to create new readout.mod 
%      and tipdown.mod files


%% update modules.txt

% read in modules.txt file as string
full_text = fileread('modules.txt');

% replace given modules with correct names
full_text = strrep(full_text,sprintf('module%d.mod',tipdown_index),'tipdown.mod');
full_text = strrep(full_text,sprintf('module%d.mod',readout_index),'readout.mod');

% rewrite modules.txt
fid = fopen('modules.txt','wt');
fprintf(fid, full_text);
fclose(fid);


%% make copies of .mod files with correct names
system(['cp module',num2str(tipdown_index),'.mod tipdown.mod']);
system(['cp module',num2str(readout_index),'.mod readout.mod']);



end

