% To use the +toppe Matlab package, include the (root) folder containing the +toppe folder in your Matlab path.
addpath ~/gitlab/toppe/

% step
% 0 ... Basic sequence
% 1 ... Add spoiler in read, phase and slice (vary spoiler?)
% 2 ... Refocus in phase
% 3 ... Vary RF phase quasi-randomly
% 4 ... Make receiver phase follow transmitter phase
step = 4;

% Define FOV and resolution
n = 128;
nz = 1;
matrix = [n n nz];
sliceThickness = 0.5;    % cm
fov  = 25.6;             % cm

% Define sequence parameters
alpha = 20;              % excitation angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x)

% set system limits
sys = toppe.systemspecs('maxSlew', 130, 'slewUnit', 'T/m/s', 'maxGrad', 25, 'gradUnit', 'mT/m');  

% Create slice selection pulse and gradient ('tipdown.mod')
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
ftype = 'ls';
tbw = 8;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(alpha, sliceThickness, tbw, dur, ncyclesspoil, ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
%toppe.plotmod(ofname);

% Create readout.mod
ofname = 'readout.mod';     % output file name
zres = 10*n/matrix(1);      % 'dummy' z resolution since we're only doing 2D here. Just needs to be larger than in-plane resolution.
toppe.utils.makegre(fov, matrix(1), zres, ... 
                    'system', sys, 'ofname', ofname, 'ncycles', ncyclesspoil); 
toppe.plotmod(ofname);

%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ii = 1;                 % counts number of module executions

toppe.write2loop('setup');
for iy = -10:n   % We'll use iy<1 for approach to steady-state

	% rf excitation 
	toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'Gamplitude', [0 0 1]');

	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
	daqphs = rfphs;
	yamp = max( ((iy-1+0.5)-n/2)/(n/2), -1);    % phase-encode amplitude scaling, range is [-1 1]
	toppe.write2loop('readout.mod', 'DAQphase', daqphs, 'view', max(iy,1), ...
		'Gamplitude', [1 yamp 0]', 'textra', 20); 

	% update rf/rec phase (RF spoiling)
	rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
	rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
end
toppe.write2loop('finish');


%% Play sequence in loop (movie) mode
nModulesPerTR = 2;
toppe.playseq(nModulesPerTR, 'nTRskip', 2);

return;

