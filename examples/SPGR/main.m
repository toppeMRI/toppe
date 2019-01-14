function main
% 3D SPGR sequence example

%% Set system limits
% Good to back off a bit on slew to avoid PNS (even if scanner can slew at higher rate).
% Otherwise, accept default values.
sys = toppe.systemspecs('maxSlew', 10, 'maxGrad', 4);  

%% Acquisition parameters
% fov and voxel size must be square (in-plane)
matrix = [240 240 10];
fov  = [24 24 10];       % cm
flip = 10;               % excitation flip angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

if matrix(1) ~= matrix(2) | fov(1) ~= fov(2)
	error('In-plane fov and voxel size must be square');
end

%% Make .mod file containing z spoiler gradient and RF excitation
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
slthick = fov(3)*0.8;       % Slab thickness (cm). A bit smaller than fov(3) to avoid aliasing.
ftype = 'min';              % minimum-phase SLR pulse (good for 3D imaging)
tbw = 8;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncyclesspoil*matrix(3), ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
%toppe.plotmod(ofname);

%% Create readout.mod
ofname = 'readout.mod';     % output file name
toppe.utils.makegre(fov(1), matrix(1), fov(3)/matrix(3), ... 
                    'system', sys, 'ofname', ofname, 'ncycles', ncyclesspoil); 
%toppe.plotmod(ofname);

%% Create scanloop.txt
rfmod = 1;           % module index, i.e., line number in modules.txt
readoutmod = 2;
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ii = 1;                 % counts number of module executions
ny = matrix(2);
nz = matrix(3);

dabon = 1;
daboff = 0;
dabmode = dabon;   % best to leave acquisition on even during disdaqs so auto-prescan gets a signal (?)
waveform = 1;
textra = 0;        % add delay at end of module (int, microseconds)

for iz = 0:nz           % We'll use iz=0 for approach to steady-state
	for iy = 1:ny

		% rf excitation block (usage of 'block' here parallels its usage in Pulseq)
		block = [];
		block.module = rfmod;
		block.rfscale = 1.0;
		block.gxscale = 0;
		block.gyscale = 0;
		block.gzscale = 1.0;
		block.rfphs = angle(exp(1i*rfphs));   % radians
		d(ii,:) = toppe.blockstruct2vec(block);
		ii = ii + 1;

		% readout
		block = [];
		block.module =  readoutmod;
		block.rfscale = 1.0;
		block.gxscale = 1.0;
		block.gyscale = ((iy-1+0.5)-ny/2)/(ny/2);    % phase-encode amplitude scaling, range is [-1 1]
		block.gzscale = ((iz-1+0.5)-nz/2)/(nz/2);    % partition-encode amplitude scaling
		block.dabslice = max(iz,0);                  % Convention: skip dabslice=0 
		block.dabecho = 0; 
		block.dabview = iy;                          % Convention: skip baseline (0) view
		block.recphs = angle(exp(1i*rfphs));         % radians
		d(ii,:) = toppe.blockstruct2vec(block);
		ii = ii + 1;

		% update rf/rec phase
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end

% write loop array to scanloop.txt
toppe.loop2txt(d);

%% Play sequence in loop (movie) mode
nModulesPerTR = 2;
toppe.playseq(nModulesPerTR, 'nTRskip', 10);

return;

