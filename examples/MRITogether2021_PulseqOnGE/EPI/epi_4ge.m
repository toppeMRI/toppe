% Create 3D EPI variable flip angle SPGR sequence for T1 mapping.
% Acquires multiple 3D EPI image volume, each with a different flip angle.

%% System hardware specs
% It is recommended that maxRF match the physical scanner limits, 
% to ensure accurate B1 scaling.
% maxGrad and maxSlew can be <= physical limit.
seq.sys = toppe.systemspecs('maxRF', 0.25, 'rfUnit', 'Gauss', ...
    'maxSlew', 9, 'slewUnit', 'Gauss/cm/ms', ... % to keep PNS low
    'maxGrad', 5, 'gradUnit', 'Gauss/cm', ...
    'timessi', 100);   % us

%% Readout parameters
% Must be chosen so the resulting echo-train length is an integer.
nz = 60;
seq.fov = [24 24 0.1*nz];         % cm
seq.matrix = [240 240 nz];
seq.nshots = 40; % number of in-plane EPI shots 
seq.Ry = 1;  % ky undersampling factor
seq.flyback = true; 
seq.rampsamp = false;
seq.res = seq.fov./seq.matrix;   % voxel size (cm)
seq.flip = [5 10 20 30]; % degrees

%% Create slab excitation module (tipdown.mod)
seq.rf.slThick = 0.50*seq.fov(3);       % cm
seq.rf.tbw = 6;              % time-bandwidth product of SLR pulse 
seq.rf.dur = 1.5;            % pulse duration (ms)
seq.rf.ftype = 'ls';  % easier to balance than minimum phase pulse ('min')


% common to all sequences in this folder (see main.m)
seq.mods.readout = 'readout.mod';    % see toppe.utils.makeepi()
seq.mods.prephaser = 'prephaser.mod';  % see toppe.utils.makeepi()

% Create readout.mod and prephaser.mod, which will be used by all sequences
[gx,gy,gz] = toppe.utils.makeepi(seq.fov, seq.matrix, ...
    seq.nshots, seq.sys, ...
    'flyback', seq.flyback, ...
    'rampsamp', seq.rampsamp, ...
    'Ry', seq.Ry, ...
    'isbalanced', true, ...
    'writefiles', true); 

% echo spacing (ms)
tsamp = seq.sys.raster*1e3;   % gradient sample time (ms)
if seq.flyback
    seq.es = tsamp * (numel(gx.echo) + numel(gx.flyback));
else
    seq.es = tsamp * numel(gx.echo);
end
if seq.nshots == seq.matrix(2)
    seq.es = 0;
end



%% bSSFP parameters
switch scandesign
    % times in ms
    % angles in degrees
    case 1  % D20, see email 10/29/2021
        seq.bssfp.alpha = [10.0 10.0 40.0 10.0 10.0 32.7 10.0 10.0 40.0 10.0 10.0 31.8 10.0 40.0 10.0 10.0 32.6 10.0 10.0 34.5];
        seq.bssfp.phi = [-171.9 -143.9 -135.9 -120.3 -93.4 -88.9 -67.1 -34.2 -27.3 -10.0 17.3 21.0 45.2 81.7 82.4 111.8 122.4 135.0 161.7 169.3];
        seq.bssfp.tr = 20*ones(size(seq.bssfp.alpha));
    case 2  % D40, see email 10/29/2021
        seq.bssfp.alpha = [10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 40.0 10.0 10.0 40.0 10.0 40.0 10.0 40.0];
        seq.bssfp.phi = [-176.4 -168.8 -159.5 -150.3 -142.1 -130.1 -124.4 -111.5 -107.6 -93.2 -90.5 -74.2 -73.6 -56.1 -54.7 -39.4 -37.2 -22.5 -18.0 -5.3 1.3 11.6 18.8 28.9 38.6 45.8 57.9 63.1 76.5 79.9 95.2 97.0 113.3 113.9 131.3 133.3 148.5 153.1 166.1 172.1];
        seq.bssfp.tr = 20*ones(size(seq.bssfp.alpha));
    case 3  % test scan
        seq.bssfp.alpha = [40 10];  % high flip first for setting receive gain
        seq.bssfp.phi = [180 180];
        seq.bssfp.tr = 20*ones(size(seq.bssfp.alpha));
end

seq.bssfp.flip = max(seq.bssfp.alpha(:));  % max flip angle

