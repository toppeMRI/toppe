function [delay90 delay180] = getsedelays(sys, tipdown, refocus, readout, te)
% Calculate delays needed to achieve a desired TE in a spin-echo train.
% Sequence is assumed to consist of tipdown-refocus-readout-refocus-readout-refocus...

% Inputs:
%  tipdown   [string]   .mod file name containing 90 degree excitation pulse
%  refocus   [string]   .mod file containing 180 refocusing pulse
%  readout   [string]   .mod file containing readout
%  te        [1 1]      desired TE (ms)
%
% Outputs:
%  delay90   [1 1]   (ms). Apply after 90 pulse
%  delay180  [1 1]   (ms). Apply before AND after 180 pulses
%

dt = 4e-3;    % ms

% get minimum TE
toppe.write2loop('setup', sys);
toppe.write2loop(tipdown, sys);
toppe.write2loop(refocus, sys);
toppe.write2loop(readout, sys);
toppe.write2loop(refocus, sys);
toppe.write2loop('finish', sys);
minsep180 = toppe.getTRtime(2, 3, sys)*1e3;     % ms. Minimum separation between 180 pulses.

% return required delay of 180 pulse (to achieve desired TE)
if minsep180 > te
	error('Minimum separation between 180 refocusing pulses exceeds requested TE');
else
	delay180 = (te-minsep180)/2;
end

% get separation of 90 pulse and first 180 pulse
rf = toppe.plotseq(1, 1, sys, 'doDisplay', false);
n1 = length(rf);                                 % length of 90 module
rf = toppe.plotseq(2, 2, sys, 'doDisplay', false);
n2 = length(rf);                                 % length of spin-echo module

rf = toppe.plotseq(1, 2, sys, 'doDisplay', false);
i90 = find(abs(rf(1:n1))==max(abs(rf(1:n1))));
i90 = i90(1);                                    % the peak can contain more than one sample
I = find(abs(rf(n1:end))==max(abs(rf(n1:end))));
i180 = I(1) + n1;                               
minsep90180 = dt*(i180-i90);
if minsep90180 > te/2
	error('Minimum separation between 90 and first 180 pulse exceeds requested TE');
else
	delay90 = te/2 - minsep90180;
end
