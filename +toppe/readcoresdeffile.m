function cores = readcorelistfile(fn)
% function cores = readModuleListFile(fn='cores.txt')
%
% Input:
%  fn           default: 'cores.txt'
%
% Output:
%  cores        cell array where each entry contains the module IDs for one row from cores.txt

if ~exist('fn', 'var')
	fn = 'cores.txt';
end

fid = fopen(fn, 'r', 'ieee-be');

s = fgets(fid);  % skip line

ncores = fscanf(fid, '%d\n', 1);

s = fgets(fid);  % skip line

for ic = 1:ncore
    % number of module instances (columns in cores.txt) for this core
    nModInstance = fscanf(fid, '%d', 1); 

    cores{ic} = fscanf(fid, '%d', nModInstance);  % [1 nModInstance]
end
 
fclose(fid);

