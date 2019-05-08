% 2D multi-shot EPI demo using TOPPE @ ISMRM 2019
%
% Jon-Fredrik Nielsen, University of Michigan 

% To use the +toppe Matlab package, include the (root) folder containing the +toppe folder in your Matlab path.
addpath ../..   % path to the 'toppe' folder

%% Define sequence parameters
N = 96;                  % image size (pixels)
fov = 25.6;              % field of view (cm)
slicethickness = 0.5;    % cm
nshots = 8;              % number of RF shots to fill (2D) k-space
alpha = 40;              % excitation angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

%% set system limits
% NB! 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' is a design choice, i.e., it can be at or below the physical system limit.
sys = toppe.systemspecs('maxSlew', 150, 'slewUnit', 'T/m/s', 'maxGrad', 50, 'gradUnit', 'mT/m');  

%% Create slice selection pulse and gradient ('tipdown.mod')
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
ftype = 'ls';               % ls = least-squares SLR design
tbw = 6;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(alpha, slicethickness, tbw, dur, ncyclesspoil, ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
%toppe.plotmod(ofname);

%% Create readout.mod
ofname = 'readout.mod';     % output file name
toppe.utils.makeepi(fov,N,nshots,'ofname',ofname, 'system', sys);
%toppe.plotmod(ofname);

%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees

[~,~,~,~,~,hdrints] = toppe.readmod('readout.mod');
raster = 4e-3;                     % msec
deltaTE = raster*hdrints(3);       % echo spacing (msec)

toppe.write2loop('setup');
for iy = -2:nshots   % We'll use iy<1 for approach to steady-state

	% rf excitation 
	ets = (max(iy,1)-1)/nshots*deltaTE;     % echo-time shift (msec)
	toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'Gamplitude', [0 0 1]', 'textra', ets);

	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
	tdelay = 960-ets;     % (ms) delay at end of waveforms (determines TR)
	toppe.write2loop('readout.mod', 'waveform', max(iy,1), 'DAQphase', rfphs, 'view', max(iy,1), ...
		'Gamplitude', [1 1 0]', 'textra', tdelay); 

	% update rf phase (RF spoiling)
	rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
	rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
end
toppe.write2loop('finish');

%% Create archive with all scan files
system('tar czf epi.tgz modules.txt scanloop.txt tipdown.mod readout.mod');

%% Play sequence in loop (movie) mode
nModulesPerTR = 2;
toppe.playseq(nModulesPerTR, 'tpause', 0.66, 'drawpause', false);

return;

