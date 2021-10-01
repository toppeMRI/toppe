function playseq(nModPerTR, system, varargin)
% Display TOPPE sequence in movie/loop mode.
%
% function playseq(nModPerTR, system, varargin)
%
% INPUTS:
%    nModPerTR          (int) number of modules per TR
%    system 			System struct
% Options:
%    nTRskip            Display only every nTRskip TRs (for speeding up loop) (default: 0)
%    loopFile           Default: 'scanloop.txt'
%    moduleListFile     Default: 'modules.txt'
%    tpause             Delay before displaying next TR (sec) (Default: 0.01)
%    drawpause          Display pauses (zeros at end of plot) (Default: 1)
%    gmax               Gauss/cm
%    rhomax             Gauss

import toppe.*
import toppe.utils.*

if nargin < 1
	error('nModPerTR must be specified');
end
if nargin < 2
    error('missing system argument');
end

%% parse inputs
% Default values 
arg.nTRskip         = 0;
arg.loopFile        = 'scanloop.txt';
arg.moduleListFile  = 'modules.txt';
arg.tpause          = 0.01; 
arg.drawpause       = 1;
arg.gmax            = 5;     % Gauss/cm. Display limit.
arg.rhomax          = 0.25;  % Gauss. Display limit.

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% Load files and display
% read scanloop
loopArr = toppe.tryread(@toppe.readloop, arg.loopFile);

% read module waveforms
modules = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);

% display
for ii = 1 : ((1+arg.nTRskip)*nModPerTR) : size(loopArr,1)
	toppe.plotseq(ii, ii+nModPerTR-1, system, 'loopArr', loopArr, 'mods', modules, ...
		'drawpause', arg.drawpause, ...
		'gmax', arg.gmax, 'rhomax', arg.rhomax );
	subplot(511); title(num2str(ii));
	drawnow; pause(arg.tpause);     % Force display update and pause
end
