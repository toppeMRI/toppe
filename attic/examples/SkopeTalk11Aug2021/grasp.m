function grasp
% Golden-angle stack of stars dynamic 3D imaging

%% Acquisition parameters
res = 0.1;    % spatial resolution (cm) (isotropic)
nz = 140; 
fov = [22 22 res*nz];   % cm
nx = 2*round(fov(1)/res/2);

flip = 5;                   % degrees
slabThick = fov(3)*0.8;      % cm

sliceOffset = 0;  % cm

scanDur = 60;                        % approximate total scan duration (sec)
TR = 5e-3;                           % approximate TR (sec)
nSpokes = 2*round(scanDur/TR/nz/2);  % number of spokes per kz encoding

%% Create modules.txt
% Entries are tab-separated.
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname	duration(us)	hasRF?	hasDAQ?\n' ...
'readout.mod	0	0	1\n' ...
'tipdown.mod	0	1	0'];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

%% set system limit struct
% NB! If passing 'sys' to writemod.m, 'maxGrad' MUST match the physical system
% limit -- since gradients are scaled relative to this.
%d 'maxSlew' can always be a design choice, i.e., it can be at or below the 
% physical system limit.
sys = toppe.systemspecs('maxSlew', 15, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm'); 

%% Create slab-selective excitation module (tipdown.mod)
% NB! Instead of calling toppe.utils.rf.makeslr, 
% the use of sigPy.mri.rf is recommended
tbw = 12;                    % time-bandwidth product
dur = 0.7;                   % ms
nCycles = 1.5 * nz;    % cycles of phase of spoiling over excited slab
ftype = 'min';
type = 'st';         % 'st' = small-tip. 'ex' = 90 degree design
[rf,gex,freqOffset] = toppe.utils.rf.makeslr(flip, slabThick, ...
    tbw, dur, nCycles, ...
    'ofname', 'tipdown.mod', ...
    'ftype', ftype, ...
    'type', type, ...
    'sliceOffset', sliceOffset, ...
    'spoilDerate', 0.5, ...
    'system', sys);

%% Create stack of stars readout module (readout.mod)
dz = res;   % z voxel size (cm)
nCycleSpoil = 0;  % balanced
oprbw = 125/2;
toppe.utils.makegre(fov(1), nx, dz, ...
    'system', sys, ...
    'ncycles', nCycleSpoil, ...
    'slewDerate', 0.6, ...   % reduce PNS
    'ofname', 'tmp.mod', ...
    'oprbw', oprbw);
[~,gx,gy,gz,desc,hdrints] = toppe.readmod('tmp.mod');
system('rm tmp.mod');
toppe.writemod('ofname', 'readout.mod', 'gx', gx, 'gy', 0*gx, 'gz', gz, ...
    'desc', desc, 'hdrints', hdrints, 'system', sys);

%% Create scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;

% golden angle (radians). From Lingala et al MRM 77:112-125 (2017)
ga = 2*pi*2/(sqrt(5)+1);

phi = 0;             % in-plane spiral rotation angle [radians]
toppe.write2loop('setup', 'version', 3);
for iview = 0:nSpokes
    for iz = 1:nz
        % z encoding gradient amplitude (normalized)
        a_gz = (iz > 0) * ((iz-1+0.5)-nz/2)/(nz/2);

		% rf excitation
		toppe.write2loop('tipdown.mod', 'RFphase', rfphs, 'RFoffset', round(freqOffset));

		% readout. Data is stored in 'slice', 'echo', and 'view' indeces
		toppe.write2loop('readout.mod', ...
            'Gamplitude', [1 1 -a_gz]', ...
            'DAQphase', rfphs, ...
            'slice', max(iz, 1), ...
            'view', max(iview, 1), ...
            'dabmode', 'on', ...
            'rot', phi);  % in-plane rotation angle (rad)

		% update rf phase (RF spoiling)
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;

	    % spiral staircase CAIPI sampling along kz
        phi = phi + 2*pi/nz;
	end

	% update rotation angle
	phi = phi + ga;
end
toppe.write2loop('finish');

TR = toppe.getTRtime(1,2);
fprintf('TR: %.2f ms\n', TR*1e3);
fprintf('Scan duration: %.1f s\n', toppe.getscantime());

%% create tar file
system('tar czf scan,grasp.tgz grasp.m tipdown.mod readout.mod modules.txt scanloop.txt');

%% display sequence
%toppe.playseq(2, 'tpause', 0.5);

% display .mod files
toppe.plotmod('all');

% display first two lines in scanloop.txt
figure; toppe.plotseq(1,2);

return;


