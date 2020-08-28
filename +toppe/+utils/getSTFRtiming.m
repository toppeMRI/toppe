function [Tfree Tg TE] = getSTFRtiming(TEro, varargin)
% function [Tfree Tg TE] = getSTFRtiming(TEro, varargin)
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
% Input options:
%  tipdown   [string]   .mod file name containing tipdown excitation pulse. Default: 'tipdown.mod'
%  readout   [string]   .mod file containing readout. Default: 'readout.mod'
%  tipup     [string]   .mod file name containing tipup excitation pulse. Default: 'tipup.mod'
%  sys       struct     System hardware specs. See toppe.systemspecs. Default: sys = toppe.systemspecs();
%
% Outputs:
%  TR        double     Sequence TR (ms)
%  Tfree     double     Time from peak of tipdown pulse to peak of tipup pulse (ms)
%  Tg        double     Time from peak of tipup pulse to peak of tipdown pulse (ms)
%  TE        double     Time to echo (ms)
%
% Example:
% >> [~,~,~,~,~,hdrints] = toppe.readmod('readout.mod');
% >> sys = toppe.systemspecs();
% >> TEro = 1e-3*(sys.toppe.start_core_daq + sys.toppe.daqdel + hdrints(3)*4);  % ms
% >> [TR, Tfree, Tg, TE] = toppe.utils.getSTFRtiming(TEro, 'tipdown', 'tipdown,1.mod', 'readout', 'readout,1.mod', 'tipup', 'tipup,1.mod', 'spoiler', 'spoiler,1.mod');


% defaults
arg.tipdown = 'tipdown.mod';
arg.readout = 'readout.mod';
arg.tipup   = 'tipup.mod';
arg.spoiler = 'spoiler.mod';
arg.system     = toppe.systemspecs();

% parse inputs
arg = toppe.utils.vararg_pair(arg, varargin);

dt = 4e-3;    % ms

% write scanloop.txt
toppe.write2loop('setup');
toppe.write2loop(arg.tipdown);
toppe.write2loop(arg.readout);
toppe.write2loop(arg.tipup);
toppe.write2loop(arg.spoiler);
toppe.write2loop('finish');

% Tfree
rf = toppe.plotseq(1, 1, 'doDisplay', false, 'system', arg.system);
n1 = length(rf);                                 % length of tipdown 
rf = toppe.plotseq(2, 3, 'doDisplay', false, 'system', arg.system);
n2 = length(rf);                                 % length of readout+tipup 
rf = toppe.plotseq(1, 3, 'doDisplay', false, 'system', arg.system);
I = find(abs(rf(1:n1))==max(abs(rf(1:n1))));
itipdown = I(1);                                    % the peak can contain more than one sample
I = find(abs(rf(n1:end))==max(abs(rf(n1:end))));
itipup = I(1) + n1;                               
Tfree = dt*(itipup-itipdown);

tipdown.dur = dt*n1;            % ms
tipdown.tpeak = dt*itipdown;    % ms

% Tg
rf = toppe.plotseq(1, 4, 'doDisplay', false, 'system', arg.system);
Tg = dt*length(rf) - Tfree;

% TE
TE = (tipdown.dur-tipdown.tpeak) + TEro;
