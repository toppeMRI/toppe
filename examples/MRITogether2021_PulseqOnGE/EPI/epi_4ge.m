% Create 3D EPI variable flip angle SPGR sequence for T1 mapping.
% Acquires multiple 3D EPI image volume, each with a different flip angle.

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
seq.fov = [24 24 6];         % cm
seq.matrix = 10*seq.fov;     % isotropic resolution
seq.nshots = 40;             % number of EPI segments
if rem(seq.matrix(2), seq.nshots) ~= 0
    error('Echo-train length must be an integer');
end
seq.Ry = 1;                  % ky undersampling factor
seq.flyback = true; 
seq.rampsamp = false;
seq.res = seq.fov./seq.matrix;   % voxel size (cm)
seq.flip = [5 10 20 30];         % degrees


%% Create slab excitation module (ex.mod)
% In TOPPE, all waveforms have 4us raster time (both RF and gradients)

% Design waveforms.
% We will use a wrapper for John Pauly's SLR toolbox for this.
seq.rf.slThick = 0.80*seq.fov(3);     % cm
seq.rf.tbw = 6;              % time-bandwidth product of SLR pulse 
seq.rf.dur = 1.5;            % pulse duration (ms)
seq.rf.ftype = 'min';        %  minimum-phase SLR design
nCyclesSpoil = 2;  % cycles of gradient spoiling along z (slice-select direction)
[rf, gz] = toppe.utils.rf.makeslr(max(seq.flip), seq.rf.slThick, ...
    seq.rf.tbw, seq.rf.dur, nCyclesSpoil, seq.sys, ...
    'ftype', seq.rf.ftype, ...
    'writeModFile', false);
rf = toppe.makeGElength(rf);  % make waveform length multiple of 4
gz = toppe.makeGElength(gz);  % make waveform length multiple of 4

% Write waveforms to a .mod file
toppe.writemod(seq.sys, 'rf', rf, 'gz', gz, ...
    'ofname', mods.ex);


%% Create EPI modules (prephaser.mod and readout.mod)

% Design waveforms using the makeepi helper function,
% that returns structs containing the various gradient segments.
[gx, gy, gz] = toppe.utils.makeepi(seq.fov, seq.matrix, seq.nshots, seq.sys, ...
    'flyback', seq.flyback, ...
    'rampsamp', seq.rampsamp);

% Create prephaser.mod
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


%% Calculate minimum TR (optional -- included here for instructional purposes)
toppe.write2loop('setup', seq.sys);
toppe.write2loop(mods.ex, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop(mods.readout, seq.sys);
toppe.write2loop(mods.prephaser, seq.sys);
toppe.write2loop('finish', seq.sys);

trmin = (seq.nshots-1)/seq.nshots*seq.es + toppe.getTRtime(1, 4, seq.sys) * 1e3;   % ms


%% Write scanloop.txt

ny = seq.matrix(2);
nz = seq.matrix(3);

toppe.write2loop('setup', seq.sys);
for iim = 1:nscans
    for ib = 1:30
        fprintf('\b');
    end
    fprintf('Writing scan %d of %d', iim, nscans);

    a_ex = seq.bssfp.alpha(iim)/seq.bssfp.flip;

    rfphs = 0;    % radians

    % set textra to achieve desired TR
    textra = max(0, seq.bssfp.tr(iim) - trmin);  % ms

    for iz = -3:nz   % iz < 1 are discarded acquisitions to reach steady state
        for iy = 1:seq.nshots
            % y/z phase-encode amplitudes, scaled to (-1,1)
            a_gy = -(iz>0)*((iy-1+0.5)-ny/2)/(ny/2);
            a_gz = -(iz>0)*((iz-1+0.5)-nz/2)/(nz/2);   

            % rf phase cycling
            rfphs = rfphs + seq.bssfp.phi(iim)/180*pi;  % rad

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
                'slice', iim, 'echo', max(iz,1), 'view', iy, ...
                'dabmode', 'on');

            % rephase, and add delay to achieve desired (and constant) TR
            toppe.write2loop(mods.prephaser, seq.sys, ...
                'Gamplitude', [-1 -a_gy -a_gz]', ...
                'textra',   seq.es - iy/seq.nshots*seq.es + textra);
      end
    end
end
fprintf('\n');
toppe.write2loop('finish', seq.sys);


%% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in toppeN.entry
% NB! The file toppeN.entry must exist in the folder from where this script is called.
toppe.preflightcheck('toppeN.entry', 'seqstamp.txt', seq.sys);

%% create tar file
system('tar czf ~/tmp/scan,exchange,bssfp.tgz seqstamp.txt modules.txt scanloop.txt *.mod getparams.m bssfp.m');

return;

%% simulate slice profile
[rf,~,~,gz] = toppe.readmod('ex.mod');
X = linspace(-seq.rf.slThick, seq.rf.slThick, 100);         % cm
dt = 4e-3;            % ms
Gamp = max(gz);       % Gauss
T1 = 500;
T2 = 50;
toppe.utils.rf.slicesim(rf,dt,T1,T2,Gamp,X);

