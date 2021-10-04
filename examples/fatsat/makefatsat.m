% Create fatsat.mod
%
% Use makeslr() function, which is a wrapper for John Pauly's toolbox,
% to create a pulse with bandwidth 500 Hz.
%
% To use this pulse in a TOPPE scan, set the 'RFoffset' parameter to -440 (for 3T) in
% the scan loop file, and turn off the gradients:
%
%   toppe.write2loop('fatsat.mod', ...
%                    'Gamplitude', [0 0 0]', ... 
%                    'RFoffset', -440);

% system hardware specs
sys = toppe.systemspecs('maxSlew', 20, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'maxRF', 0.25);  % Gauss

% design SLR pulse and write to 'fatsat.mod'
flip = 90;
tbw = 1.5;       % time-bandwidth product
dur = 3;         % pulse duration (msec)
toppe.utils.rf.makeslr(flip, 1e5, tbw, dur, 1e-8, sys, ...
    'ftype', 'ls', ...
    'type', 'ex', ...
    'ofname', 'fatsat.mod');

% plot
toppe.plotmod('fatsat.mod');
