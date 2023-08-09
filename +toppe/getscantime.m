function dur = getscantime(sysGE, varargin)
% function dur = getscantime(system, varargin)
%
% Get total scan time for a TOPPE scan.
% Like the interpreter on the scanner, this function
% load the sequence files present in the current folder.
%
% Input:
%   sysGE              struct specifying hardware system info, see systemspecs.m
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

%% Call plotseq with 'doTimeOnly' argument to get the sequence duration quickly
[rf,gx,gy,gz,T] = toppe.plotseq(sysGE, 'timeRange', [0 inf], 'doTimeOnly', true);
dur = T(2);
%fprintf('Total scan time: %dm %0.1f s\n', floor(dur/60), dur - 60*floor(dur/60) );
