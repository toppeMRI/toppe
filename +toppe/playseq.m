function playseq(nModPerTR,varargin)
% Display TOPPE sequence in movie/loop mode.
%
% function playseq(nModPerTR,varargin)
%
% INPUTS:
%    nModPerTR      number of modules per TR
% Options:
%    nTRskip            Display only every nTRskip TRs (for speeding up loop) (default: 0)
%    loopFile           Default: 'scanloop.txt'
%    moduleListFile    Default: 'modules.txt'
%    system 				system struct. Default calls systemspecs.m.
%    tpause             Delay before displaying next TR (sec) (default: 0)

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
% (c) 2016-18 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: playseq.m,v 1.4 2018/10/25 12:40:29 jfnielse Exp $

import toppe.*
import toppe.utils.*

if nargin < 1
	error('nModPerTR must be specified');
end
if ~strcmp(class(nModPerTR), 'double') | nModPerTR < 1
	error('First argument (nModPerTR) must be a number.');
end

%% parse inputs
% Default values 
arg.nTRskip         = 0;
arg.loopFile        = 'scanloop.txt';
arg.moduleListFile = 'modules.txt';
arg.system          = toppe.systemspecs();
arg.tpause          = 0.01; 

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% Load files and display
% read scanloop
loopArr = toppe.utils.tryread(@toppe.readloop, arg.loopFile);

% read module waveforms
modules = toppe.utils.tryread(@toppe.readmodulelistfile, arg.moduleListFile);

% display
for ii = 1 : ((1+arg.nTRskip)*nModPerTR) : size(loopArr,1)
	toppe.plotseq(ii, ii+nModPerTR-1, 'loopArr', loopArr, 'mods', modules, 'system', arg.system);
	subplot(511); title(num2str(ii));
	pause(arg.tpause);     % to allow display to refresh
end
