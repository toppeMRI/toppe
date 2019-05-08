% To use the +toppe Matlab package, include the (root) folder containing the +toppe folder in your Matlab path.
addpath ../..   % path to the 'toppe' folder
addpath ./lib   % files needed for SLR design 

% step
% 0 ... Basic sequence
% 1 ... Add RF spoiling
% 2 ... Make receiver phase follow transmitter phase
step = 2;

% Define FOV and resolution
n = 128;
nz = 1;
matrix = [n n nz];
sliceThickness = 0.5;    % cm
fov = 25.6;              % cm

% Define sequence parameters
alpha = 20;              % excitation angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

% set system limits
% NB! 'maxGrad' MUST match the physical system limit -- since gradients are scaled relative to this.
% 'maxSlew' can be a design choice, i.e., it can be at or below the physical system limit.
sys = toppe.systemspecs('maxSlew', 150, 'slewUnit', 'T/m/s', 'maxGrad', 50, 'gradUnit', 'mT/m');  

% Create slice selection pulse and gradient ('tipdown.mod')
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
ftype = 'ls';
tbw = 6;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(alpha, sliceThickness, tbw, dur, ncyclesspoil, ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
%toppe.plotmod(ofname);

% Create readout.mod
ofname = 'readout.mod';     % output file name
zres = 10*n/matrix(1);      % 'dummy' z resolution since we're only doing 2D here. Just needs to be larger than in-plane resolution.
toppe.utils.makegre(fov, matrix(1), zres, ... 
                    'system', sys, 'ofname', ofname, 'ncycles', ncyclesspoil); 
%toppe.plotmod(ofname);

%% Create scanloop.txt
rfphs = 0;              % radians
daqphs = 0;
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees

toppe.write2loop('setup');
for iy = -10:n   % We'll use iy<1 for approach to steady-state

	% rf excitation 
	toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'Gamplitude', [0 0 1]');

	% readout. Data is stored in 'slice', 'echo', and 'view' indeces.
	if step > 1
		daqphs = rfphs;
	end
	yamp = max( ((iy-1+0.5)-n/2)/(n/2), -1);    % phase-encode amplitude scaling, range is [-1 1]
	toppe.write2loop('readout.mod', 'DAQphase', daqphs, 'view', max(iy,1), ...
		'Gamplitude', [1 yamp 0]', 'textra', 30); 

	% update rf phase (RF spoiling)
	if step > 0
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end
toppe.write2loop('finish');

%% Create archive with all scan files
system('tar czf gre.tgz modules.txt scanloop.txt tipdown.mod readout.mod');


%% Play sequence in loop (movie) mode
nModulesPerTR = 2;
toppe.playseq(nModulesPerTR, 'nTRskip', 2);

return;

