function [d, hdr] = readloop(loopfile)
% read scanloop.txt (toppe.e scan loop definition, used with modules.txt).
%
% function d = readloop(loopfile)
%
% Input:
%    loopfile       default: 'scanloop.txt'

% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
% TOPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU Library General Public License as published by
% the Free Software Foundation version 2.0 of the License.
%
% TOPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public License
% along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
% 
% (c) 2016 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

if ~exist('loopfile', 'var')
	loopfile = 'scanloop.txt';
end

% read header
fid = fopen(loopfile, 'r', 'ieee-be');
s = fgets(fid);  % read a whole line (here we just skip it)
hdr.nt       = fscanf(fid, '%d\t', 1);
hdr.maxslice = fscanf(fid, '%d\t', 1);
hdr.maxecho  = fscanf(fid, '%d\t', 1);
hdr.maxview  = fscanf(fid, '%d\t', 1);
hdr.scandur  = fscanf(fid, '%d\t', 1);
hdr.version  = fscanf(fid, '%d\n', 1);
fclose(fid);

% read data (tab-separated, starting at row 4)
d = dlmread(loopfile,'\t',3,0);

return;
