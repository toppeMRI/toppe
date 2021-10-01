function dur = getTRtime(LineStart,LineEnd,system,varargin)
% Get scan time to run lines LineStart through LineEnd of a scanloop.
% Useful for calculating the TR of just a portion of your scan
%
% function dur = getTRtime(LineStart,LineEnd,system,varargin)
%
%
% Input options:
%   loopFile           Default: 'scanloop.txt'
%   moduleListFile     Default: 'modules.txt'
%   system             struct specifying hardware system info, see systemspecs.m
%                      Default: arg.system = systemspecs();
%   In addition, the .mod files listed in 'moduleListFile' must be present in the current (working) folder.
%
% Output:
%    dur          sec
%
% Examples:
% >> toppe.getTRtime(100,102); %Time to run lines 100 through 102

import toppe.*
import toppe.utils.*

%% parse inputs
% Default values 
arg.loopFile       = 'scanloop.txt';
arg.moduleListFile = 'modules.txt';

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% read scan files
% read scanloop
loopArr = toppe.tryread(@toppe.readloop, arg.loopFile);

% read module content
mods = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);

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
for ii = LineStart:LineEnd
	rho = toppe.plotseq(ii, ii, system, 'loopArr', loopArr, 'mods', mods, 'doDisplay', false);
	dur = dur + size(rho,1)*dt;
end

