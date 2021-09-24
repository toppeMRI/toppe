function dur = getscantime(system, varargin)
% Get total scan time for a TOPPE scan.
%
% function dur = getscantime(system, varargin)
%
% Input:
%   system             struct specifying hardware system info, see systemspecs.m
%
% Input options:
%   loopFile           Default: 'scanloop.txt'
%   moduleListFile     Default: 'modules.txt'
%   loopArr            Directly specify contents of loopFile   
%   mods               Directly specify contents of module list
%
%   In addition, the .mod files listed in 'moduleListFile' must be present in the current (working) folder.
%
% Output:
%    dur          sec
%
% Examples:
% >> toppe.getscantime();

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
arg.loopArr        = [];
arg.mods           = [];

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% read scan files
if isempty(arg.loopArr)
    loopArr = toppe.tryread(@toppe.readloop, arg.loopFile); % read scanloop
else
    loopArr = arg.loopArr;
end

% read module content
if isempty(arg.mods)
   mods = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);
else
   mods = arg.mods;
end

%% Call plotseq with 'doTimeOnly' argument to get the length of rho quickly
% Do it in groups of rows (in scanloop.txt) to reduce memory requirements
dt = 4e-6;    % duration of one gradient/rf sample (sec)
dur = 0;
nl = size(loopArr,1);   % number of rows in scanloop.txt (startseq() calls)
nlPerGroup = 1e3;
if nl > 2*nlPerGroup
	for ii = 1:nlPerGroup:(nl-nlPerGroup)   %size(loopArr,1)
		rho = toppe.plotseq(ii, ii+nlPerGroup-1, system, 'loopArr', loopArr, 'mods', mods, 'doTimeOnly', true);
		dur = dur + size(rho,1)*dt;
	end
	nlLeft = max(0, nl - (ii+nlPerGroup-1));
	rho = toppe.plotseq(ii+1, ii+nlLeft-1, system, 'loopArr', loopArr, 'mods', mods, 'doTimeOnly', true);
	dur = dur + size(rho,1)*dt;
else
	rho = toppe.plotseq(1, size(loopArr,1), system, 'loopArr', loopArr, 'mods', mods, 'doTimeOnly', true);
	dur = dur + size(rho,1)*dt;
end
%dur = round(dur);
%fprintf('Total scan time: %dm %ds\n', floor(dur/60), dur - 60*floor(dur/60) );
