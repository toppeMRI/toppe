function modArr = readModuleListFile(modulesListFile)
% Load all .mod files in modulesListFile into cell array
%
% function modArr = readModuleListFile(modulesListFile)
%
% Input:
%  modulesListFile     default: 'modules.txt'
%
% Output:
%  modArr       cell array of structs s with members
%                   s.fname
%                   s.dur
%                   s.hasRF
%                   s.hasDAQ
%                   s.rf/gx/gy/gz/desc/ as returned by readmod.m
%
% $Id: readmodulelistfile.m,v 1.2 2018/10/24 10:44:50 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/readmodulelistfile.m,v $

if ~exist('modulesListFile', 'var')
	modulesListFile = 'modules.txt';
end

% get waveforms
fid = fopen(modulesListFile, 'r', 'ieee-be');
s = fgets(fid);  % skip line
ncores = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
for ic = 1:ncores
	modArr{ic}.fname = fscanf(fid, '%s ', 1);
	modArr{ic}.dur = fscanf(fid, '%d ', 1);
	modArr{ic}.hasRF = fscanf(fid, '%d ', 1);
	modArr{ic}.hasDAQ = fscanf(fid, '%d\n', 1);
	[modArr{ic}.rf, modArr{ic}.gx, modArr{ic}.gy, modArr{ic}.gz, ...
    desc, modArr{ic}.paramsint16, modArr{ic}.paramsfloat] = toppe.readmod(modArr{ic}.fname,false);
	modArr{ic}.wavdur = numel(modArr{ic}.gx(:,1))*4;   % waveform duration [us]
end
fclose(fid);

