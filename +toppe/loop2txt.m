function loop2txt(d)
% Write scanloop array to scanloop.txt.
%
% function loop2txt(d)
%
% Input:
%  d      [nt 16] array. See 'writeloop.m' functions in the examples.

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
% $Id: loop2txt.m,v 1.2 2018/10/24 10:44:49 jfnielse Exp $

nt = size(d,1);              % number of startseq() calls   
maxslice = max(d(:,7));
maxecho = max(d(:,8));
maxview = max(d(:,9));
fname = 'scanloop.txt';
fid = fopen(fname, 'w', 'ieee-be');
fprintf(fid, 'nt\tmaxslice\tmaxecho\tmaxview\n');
fprintf(fid, '%d\t%d\t%d\t%d\n', nt, maxslice, maxecho, maxview);
fprintf(fid, 'Core ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview dabon phi rfphase recphase \n');
fclose(fid);
dlmwrite(fname, d, '-append', 'delimiter', '\t', 'precision', 8);  % precision=8 needed to avoid large numbers written in scientific notation

return;
