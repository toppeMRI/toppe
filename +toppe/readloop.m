function d = readloop(loopfile)
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
%
% $Id: readloop.m,v 1.2 2018/10/24 11:05:10 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/readloop.m,v $

if ~exist('loopfile', 'var')
	loopfile = 'scanloop.txt';
end

%NL = 11;   % toppe2
%NL = 10;   % toppe.e and toppe_so.e
%NL = 13;   % toppe3, toppe4
%NL = 14;   % toppe5, toppe7e
%NL = 15;   % toppe8*

% read header
fid = fopen(loopfile, 'r', 'ieee-be');
s = fgets(fid);  % read a whole line
nt       = fscanf(fid, '%d\t', 1);
maxslice = fscanf(fid, '%d\t', 1);
maxecho  = fscanf(fid, '%d\t', 1);
maxview  = fscanf(fid, '%d\n', 1);
%rhrecon = fscanf(fid, '%d\n', 1);
fclose(fid);

% read data (tab-separated, starting at row 4)
d = dlmread(loopfile,'\t',3,0);

return;
