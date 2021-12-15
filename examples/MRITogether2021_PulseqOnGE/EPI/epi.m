% Create 3D EPI variable flip angle SPGR sequence for T1 mapping.
% Acquires multiple 3D EPI image volume, each with a different flip angle.
% 
% First create TOPPE files for execution on GE scanners.
% Then convert directly to a Pulseq file for execution on Siemens scanners.
%
% Repository: github/toppeMRI/toppe/examples/MRITogether2021_PulseqOnGE/EPI/


%% Write modules.txt
% We do this first since it is used later for timing calculations.
moduleListFile = 'modules.txt';

% Names of the three .mod files needed for this pulse sequence
mods.ex = 'ex.mod';                 % excitation pulse
mods.prephaser = 'prephaser.mod';   % x/y/z prephasing gradient. See toppe.utils.makeepi().
mods.readout = 'readout.mod';       % EPI readout train. See toppe.utils.makeepi().

fid = fopen(moduleListFile, 'wt');
fprintf(fid, 'Total number of unique cores\n');
fprintf(fid, '%d\n', length(fieldnames(mods)));
fprintf(fid, 'fname  duration(us)    hasRF?  hasDAQ?\n');
fprintf(fid, '%s\t0\t1\t0\n', mods.ex);
fprintf(fid, '%s\t0\t0\t0\n', mods.prephaser);
fprintf(fid, '%s\t0\t0\t1\n', mods.readout);
fclose(fid);


%% Set system hardware specs
% It is recommended that maxRF match the physical scanner limits, 
% to ensure accurate B1 scaling.
% maxGrad and maxSlew can be <= physical limit.
seq.sys = toppe.systemspecs('maxRF', 0.25, 'rfUnit', 'Gauss', ...
    'maxSlew', 9, 'slewUnit', 'Gauss/cm/ms', ... % to keep PNS low
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'timessi', 100);   % us


%% Readout parameters
seq.fov = [20 20 2];         % cm
seq.matrix = [64 64 10]; 
seq.nshots = 8;             % number of EPI segments
if rem(seq.matrix(2), seq.nshots) ~= 0
    error('Echo-train length must be an integer');
end
seq.Ry = 1;                  % ky undersampling factor
seq.flyback = true; 
seq.rampsamp = false;
seq.res = seq.fov./seq.matrix;   % voxel size (cm)
seq.flip = [5:5:40];         % degrees


%% Create slab excitation module (ex.mod)
% In TOPPE, all waveforms have 4us raster time (both RF and gradients)

% Design waveforms.
if false
% Use a wrapper for John Pauly's SLR toolbox for this.
    seq.rf.slThick = 0.80*seq.fov(3);     % cm
    seq.rf.tbw = 6;              % time-bandwidth product of SLR pulse 
    seq.rf.dur = 1.5;            % pulse duration (ms)
    seq.rf.ftype = 'min';        %  minimum-phase SLR design
    nCyclesSpoil = 2*seq.matrix(3);  % cycles of gradient spoiling across slab
    [rf, gz] = toppe.utils.rf.makeslr(max(seq.flip), seq.rf.slThick, ...
        seq.rf.tbw, seq.rf.dur, nCyclesSpoil, seq.sys, ...
        'ftype', seq.rf.ftype, ...
        'writeModFile', false);

save pulse rf gz
else
    load pulse
end

rf = toppe.makeGElength(rf);  % make waveform length multiple of 4
gz = toppe.makeGElength(gz);  % make waveform length multiple of 4

% Write waveforms to a .mod file
toppe.writemod(seq.sys, 'rf', rf, 'gz', gz, ...
    'ofname', mods.ex);


%% Create EPI modules (prephaser.mod and readout.mod)

% Design waveforms using the 'makeepi.m' helper function,
% that returns structs containing the various gradient segments.
[gx, gy, gz] = toppe.utils.makeepi(seq.fov, seq.matrix, seq.nshots, seq.sys, ...
    'flyback', seq.flyback, ...
    'rampsamp', seq.rampsamp);

% Create prephaser.mod
% The prephaser is scaled dynamically (see below), so it needs its own
% .mod file separate from the echo train.
gx.pre = toppe.makeGElength(gx.pre(:));
gy.pre = toppe.makeGElength(gy.pre(:));
gz.pre = toppe.makeGElength(gz.pre(:));
toppe.writemod(seq.sys, 'gx', gx.pre, 'gy', gy.pre, 'gz', gz.pre, ...
    'ofname', 'prephaser.mod', ...
    'desc', 'EPI prephasing gradients (move to corner of kspace)');

% Create readout.mod (EPI echo train).
% Also store a few values in header for later use (optional).
% Note that data is acquired during the entire waveform, with 4us sample/dwell time.
gx.et = toppe.makeGElength(gx.et);
gy.et = toppe.makeGElength(gy.et);
hdrints = [seq.matrix(1) seq.matrix(2) seq.nshots seq.Ry length(gx.echo)];
if seq.flyback
    hdrints = [hdrints length(gx.flyback)];
end
if ~seq.rampsamp
    hdrints = [hdrints length(gx.ramp)];
end
hdrfloats = [seq.fov(1) seq.fov(2)];

toppe.writemod(seq.sys, 'gx', gx.et, 'gy', gy.et, ...
    'ofname', 'readout.mod', ...
    'desc', 'EPI readout train', ...
    'hdrfloats', hdrfloats, ...
    'hdrints', hdrints);


%% Get echo spacing (ms)
% Needed for echo-time shift calculation
tsamp = seq.sys.raster*1e3;   % gradient sample time (ms)
if seq.flyback
    seq.es = tsamp * (numel(gx.echo) + numel(gx.flyback));
else
    seq.es = tsamp * numel(gx.echo);
end
if seq.nshots == seq.matrix(2)
    seq.es = 0;
end



%% Write scanloop.txt 
% (This will overwrite the temporary file we just created)
ny = seq.matrix(2);
nz = seq.matrix(3);

nscans = length(seq.flip);

% RF spoiling parameters
rfphs = 0;    % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees

toppe.write2loop('setup', seq.sys, 'version', 4);  % Initialize the file 
for iim = 1:nscans
    for ib = 1:50
        fprintf('\b');
    end
    fprintf('Writing scan %d of %d (flip angle = %d)', iim, nscans, seq.flip(iim));

    a_ex = seq.flip(iim)/max(seq.flip);  % scale excitation pulse

    for iz = -3:nz   % iz < 1 are discarded acquisitions to reach steady state
        for iy = 1:seq.nshots
            % y/z phase-encode amplitudes, scaled to (-1,1)
            a_gy = -(iz>0)*((iy-1+0.5)-ny/2)/(ny/2);
            a_gz = -(iz>0)*((iz-1+0.5)-nz/2)/(nz/2);   

            % rf excitation, and echo time shift
            toppe.write2loop(mods.ex, seq.sys, ...
                'RFphase', rfphs, ...
                'RFamplitude', a_ex, ...
                'textra', (iy-1)/seq.nshots*seq.es);

            % prephase (move to corner of kspace)
            toppe.write2loop(mods.prephaser, seq.sys, ...
                'Gamplitude', [1 a_gy a_gz]');

            % readout. Data is stored in 'slice', 'echo', and 'view' indeces.
            toppe.write2loop(mods.readout, seq.sys, ...
                'DAQphase', rfphs, ...
                'slice', max(iz, 1), 'echo', iim, 'view', iy, ...
                'dabmode', 'on');

            % rephase, and add delay to achieve constant TR
            toppe.write2loop(mods.prephaser, seq.sys, ...
                'Gamplitude', [-1 -a_gy -a_gz]', ...
                'textra', seq.es - iy/seq.nshots*seq.es);

            % update rf/receive phase
            rfphs = rfphs + (rf_spoil_seed/180*pi)*rf_spoil_seed_cnt ;  % radians
            rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        end
    end
end
fprintf('\n');
toppe.write2loop('finish', seq.sys);


%% Create 'sequence stamp' file for TOPPE, that the interpreter uses
%% for patient and hardware safety calculations.
% This file is listed in the 5th row in toppeN.entry
% NB! The file toppeN.entry must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', seq.sys);


%% create tar file (optional)
system('tar cf scan,epi.tar seqstamp.txt modules.txt scanloop.txt *.mod epi.m toppeN.entry README.md');


%% Create Pulseq file
addpath ~/github/toppeMRI/PulseGEq/  % +pulsegeq package
addpath ~/github/pulseq/matlab/      % +mr package

siemensLims = mr.opts('MaxGrad', seq.sys.maxGrad*10, 'GradUnit', 'mT/m', ...
    'MaxSlew', seq.sys.maxSlew*10, 'SlewUnit', 'T/m/s', ...
    'rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  

% Here we pass an empty tar file argument 
% since the TOPPE scan files already exist in the local path.
pulsegeq.ge2seq([], 'seqFile', 'epi.seq', ...
    'system', siemensLims, ...
    'systemGE', seq.sys);

%% Plot both sequences
iStart = 4*40+1; iStop = iStart + 8;  % row indeces in scanloop.txt
toppe.plotseq(iStart, iStop, seq.sys);

pulseqObj = mr.Sequence(siemensLims);
pulseqObj.read('epi.seq');
pulseqObj.plot('timeRange', [3.995 3.995+57e-3]);

return;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 Optional steps/tasks                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% simulate and plot slice profile
[rf,~,~,gz] = toppe.readmod(mods.ex);
Z = linspace(-seq.rf.slThick, seq.rf.slThick, 100);    % spatial positions (cm)
dt = 4e-3;            % ms
Gamp = max(gz);       % Gauss
T1 = 500;
T2 = 50;
m0 = [0 0 1]; % initial magnetization
m = toppe.utils.rf.slicesim(m0, rf, gz, dt, Z, T1, T2);


%% Calculate minimum TR 
% First create a short scanloop.txt file that only contains the first TR
toppe.write2loop('setup', seq.sys, 'version', 4);
toppe.write2loop(mods.ex, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop(mods.readout, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop('finish', seq.sys);

% Use the 'getTRtime' helper function to get the minimum TR,
% accounting for the time needed for echo-shifting
trmin = (seq.nshots-1)/seq.nshots*seq.es + toppe.getTRtime(1, 4, seq.sys) * 1e3   % ms

