function dur = getTRtime(LineStart,LineEnd,varargin)
% Get scan time to run lines LineStart through LineEnd of a scanloop.
% Useful for calculating the TR of just a portion of your scan
%
% function dur = getTRtime(LineStart,LineEnd,varargin)
%
%
% Input options:
%   loopFile       Default: 'scanloop.txt'
%   moduleListFile        Default: 'modules.txt'
%   system             struct specifying hardware system info, see systemspecs.m
%                      Default: arg.system = systemspecs();
%   In addition, the .mod files listed in 'moduleListFile' must be present in the current (working) folder.
%
% Output:
%    dur          sec
%
% Examples:
% >> toppe.getTRtime(100,102); %Time to run lines 100 through 102

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
% (c) 2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: getscantime.m,v 1.4 2018/10/25 12:40:29 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/getscantime.m,v $

import toppe.*
import toppe.utils.*

%% parse inputs
% Default values 
arg.loopFile       = 'scanloop.txt';
arg.moduleListFile = 'modules.txt';
arg.system         = toppe.systemspecs();

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% read scan files
% read scanloop
loopArr = toppe.utils.tryread(@toppe.readloop, arg.loopFile);

% read module content
mods = toppe.utils.tryread(@toppe.readmodulelistfile, arg.moduleListFile);

%% Check inputs
if LineStart > LineEnd
    error('Starting line must be smaller than ending line');
end
if LineEnd > size(loopArr,1)
    error('Ending line (%d) is greater than scanloop (%d).\n',LineEnd,size(loopArr,1));
end
if floor(LineStart) ~= LineStart || floor(LineEnd) ~= LineEnd
    error('Line values must be integers');
end

%% loop through scan, and tally scan duration
dt = 4e-6;    % duration of one gradient/rf sample (sec)
dur = 0;
fprintf('Calculation duration between lines %d and lines %d...',LineStart,LineEnd);
for ii = LineStart:LineEnd
	rho = toppe.plotseq(ii, ii, 'loopArr', loopArr, 'mods', mods, 'doDisplay', false, 'system', arg.system);
	dur = dur + size(rho,1)*dt;
end
fprintf(' done\n');

