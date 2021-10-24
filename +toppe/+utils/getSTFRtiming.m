function [Tfree Tg TE] = getSTFRtiming(TEro, sys, varargin)
% function [Tfree Tg TE] = getSTFRtiming(TEro, sys, varargin)
%
% Calculate TR/Tfree/Tg/TE for an STFR sequence. 
% Sequence is assumed to consist of tipdown-readout-tipup-spoiler modules.
%
% The file 'modules.txt' must exist in the current folder, and must contain the .mod files 
% passed to this function.
%
% NB! Overwrites scanloop.txt in current folder.
%
% Inputs:
%  TEro      double     Time from beginning of readout to echo (ms)
%  sys       struct     System hardware specs. See toppe.systemspecs. Default: sys = toppe.systemspecs();
% Input options:
%  tipdown   [string]   .mod file name containing tipdown excitation pulse. Default: 'tipdown.mod'
%  readout   [string]   .mod file containing readout. Default: 'readout.mod'
%  tipup     [string]   .mod file name containing tipup excitation pulse. Default: 'tipup.mod'
%  version   int        TOPPE version. Default: 4
%
% Outputs:
%  TR        double     Sequence TR (ms)
%  Tfree     double     Time from peak of tipdown pulse to peak of tipup pulse (ms)
%  Tg        double     Time from peak of tipup pulse to peak of tipdown pulse (ms)
%  TE        double     Time to echo (ms)

% defaults
arg.tipdown = 'tipdown.mod';
arg.readout = 'readout.mod';
arg.tipup   = 'tipup.mod';
arg.spoiler = 'spoiler.mod';
arg.version = 4;

% parse inputs
arg = toppe.utils.vararg_pair(arg, varargin);

dt = 4e-3;    % ms

% write scanloop.txt
toppe.write2loop('setup', sys, 'version', arg.version);
toppe.write2loop(arg.tipdown, sys);
toppe.write2loop(arg.readout, sys);
toppe.write2loop(arg.tipup, sys);
toppe.write2loop(arg.spoiler, sys);
toppe.write2loop('finish', sys);

% Tfree
rf = toppe.plotseq(1, 1, sys, 'doDisplay', false);
n1 = length(rf);                                 % length of tipdown 
rf = toppe.plotseq(2, 3, sys, 'doDisplay', false);
n2 = length(rf);                                 % length of readout+tipup 
rf = toppe.plotseq(1, 3, sys, 'doDisplay', false);
I = find(abs(rf(1:n1))==max(abs(rf(1:n1))));
itipdown = I(1);                                    % the peak can contain more than one sample
I = find(abs(rf(n1:end))==max(abs(rf(n1:end))));
itipup = I(1) + n1;                               
Tfree = dt*(itipup-itipdown);

tipdown.dur = dt*n1;            % ms
tipdown.tpeak = dt*itipdown;    % ms

% Tg
rf = toppe.plotseq(1, 4, sys, 'doDisplay', false);
Tg = dt*length(rf) - Tfree;

% TE
TE = (tipdown.dur-tipdown.tpeak) + TEro;
