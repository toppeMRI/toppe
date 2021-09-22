function modArr = readmodulelistfile(modulesListFile)
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

if ~exist('modulesListFile', 'var')
	modulesListFile = 'modules.txt';
end

% get waveforms
fid = fopen(modulesListFile, 'r', 'ieee-be');
s = fgets(fid);  % skip line
ncores = fscanf(fid, '%d\n', 1);
s = fgets(fid);  % skip line
for ic = 1:ncores
	modArr{ic}.fname = fscanf(fid, '%s\t', 1);
	modArr{ic}.dur = fscanf(fid, '%d\t', 1);
    if modArr{ic}.dur < 0
        error('Module duration must be >= 0');
    end
	modArr{ic}.hasRF = fscanf(fid, '%d\t', 1);
    if modArr{ic}.hasRF ~= 0 & modArr{ic}.hasRF ~= 1
        error('hasRF must be 0 or 1');
    end
	modArr{ic}.hasDAQ = fscanf(fid, '%d\n', 1);
    if modArr{ic}.hasDAQ ~= 0 & modArr{ic}.hasDAQ ~= 1
        error('hasDAQ must be 0 or 1');
    end
    if modArr{ic}.hasRF == 1 & modArr{ic}.hasDAQ == 1
        error('Module must be either an RF or DAQ module (or neither)');
    end
	[modArr{ic}.rf, modArr{ic}.gx, modArr{ic}.gy, modArr{ic}.gz, ...
    desc, modArr{ic}.paramsint16, modArr{ic}.paramsfloat, hdr] = toppe.readmod(modArr{ic}.fname,false);
    modArr{ic}.npre = hdr.npre;
    modArr{ic}.rfres = hdr.rfres;
	modArr{ic}.wavdur = numel(modArr{ic}.gx(:,1))*4;   % waveform duration [us]
end
fclose(fid);

