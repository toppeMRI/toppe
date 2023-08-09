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

% System hardware specs.
% Reduce slew to reduce noise and PNS
sys = toppe.systemspecs('maxSlew', 5, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'maxRF', 0.25);  % Gauss

% Design SLR pulse.
% As always in TOPPE, raster time is 4us.
flip = 90;
tbw = 1.5;       % time-bandwidth product
dur = 3;         % pulse duration (msec)
rf = toppe.utils.rf.makeslr(flip, 1e5, tbw, dur, 1e-8, sys, ...
    'ftype', 'ls', ...
    'type', 'ex', ...
    'writeModFile', 'false');
    %'ofname', 'fatsat.mod');

% Design spoiler gradient
nCycleSpoil = 4;  % number of spoiling cycles across slThick
slThick = 0.5;    % cm
gsp = toppe.utils.makecrusher(nCycleSpoil, slThick, sys, 0, sys.maxSlew, sys.maxGrad);

% Put together and write to 'fatsat.mod'
g = [0*rf; gsp];
rf = [rf; 0*gsp];
g = toppe.makeGElength(g);
rf = toppe.makeGElength(rf);
toppe.writemod(sys, 'rf', rf, 'gz', g, 'ofname', 'fatsat.mod');

% Plot
toppe.plotmod('fatsat.mod', 'gradcoil', sys.gradient);
