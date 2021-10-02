function sys = gre
% 3D RF-spoiled gradient-echo sequence example
%
% This script generates the following files:
%  modules.txt
%  scanloop.txt
%  tipdown.mod
%  readout.mod
%  seqstamp.txt
%
% These must be copied to the scanner host, in the folder specified on the first line in toppe0.entry.
% toppe0.entry (provided in this folder) is the entry point for the TOPPE interpreter,
% and must be placed in /usr/g/research/pulseq/ on the scanner host.

% Set hardware limits (for design and detailed timing calculations)
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
sys = toppe.systemspecs('maxSlew', 12.3, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm');

% Acquisition parameters
matrix = [120 120 60];
fov  = [24 24 24];       % cm
flip = 5;                % excitation flip angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

% Since we are using the helper function 'makegre' below,
% the in-plane FOV and matrix size must be square.
if matrix(1) ~= matrix(2) | fov(1) ~= fov(2)
    error('This example requires that in-plane FOV and matrix be square.');
end

% Non-selective excitation
nhard = 20; % number of waveform samples in hard pulse
nChop = [0 0]; % discarded samples before + after RF waveform. For advanced timing control.
rf = [(flip/360) / (sys.gamma * nhard * sys.raster) * ones(nhard,1)];
rf = [zeros(nChop(1)+2,1); rf; zeros(nChop(2)+2,1)];  % TOPPE wants waveforms to start and end with 0
rf = toppe.makeGElength(rf); % force number of samples to be multiple of 4
toppe.writemod(sys, ...
    'ofname', 'tipdown.mod', ...
    'nChop', nChop, ...
    'rf', rf);

% Create readout waveforms and write to readout.mod.
% Here we use the helper function 'makegre' to do that, but
% that's not a requirement.
toppe.utils.makegre(fov(1), matrix(1), fov(3)/matrix(3), sys, ... 
    'nChop', nChop, ...
    'ofname', 'readout.mod', ...
    'ncycles', ncyclesspoil); 

% Display .mod files.
%
%toppe.plotmod('all');

% Write modules.txt
modFileText = ['' ...
'Total number of unique cores\n' ...
'2\n' ...
'fname  duration(us)    hasRF?  hasDAQ?\n' ...
'tipdown.mod\t0\t1\t0\n' ...     % Entries are tab-separated   
'readout.mod\t0\t0\t1\n' ];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% Write scanloop.txt
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ny = matrix(2);
nz = matrix(3);

toppe.write2loop('setup', sys, 'version', 4);  % initialize file

for iz = 0:nz     % We'll use iz=0 for approach to steady-state
    for iy = 1:ny
        a_gy = -((iy-1+0.5)-ny/2)/(ny/2);  % negative so it starts at -kymax/-kzmax (convention)
        a_gz = -((max(iz,1)-1+0.5)-nz/2)/(nz/2);

        toppe.write2loop('tipdown.mod', sys, ...
            'RFphase', rfphs);

        toppe.write2loop('readout.mod', sys, ...
            'Gamplitude', [1.0 a_gy a_gz]', ...
            'DAQphase', rfphs, ...
            'textra', 0, ...  % slow down scan a bit (ms)
            'slice', max(iz,1), 'view', iy);

        % update rf/rec phase
        rfphs = rfphs + (rf_spoil_seed/180*pi)*rf_spoil_seed_cnt ;  % radians
        rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
    end
end

toppe.write2loop('finish', sys);  % finalize file

% Create 'sequence stamp' file for TOPPE.
% This file is listed in line 6 of toppe0.entry
toppe.preflightcheck('toppe0.entry', 'seqstamp.txt', sys);

% Write files to tar archive (for convenience only).
system('tar cf gre.tar toppe0.entry modules.txt scanloop.txt *.mod seqstamp.txt');

% Play sequence in loop (movie) mode
nModulesPerTR = 2;
toppe.playseq(nModulesPerTR, sys, ...
    'moduleListFile', 'modules.txt', ...
    'loopFile', 'scanloop.txt', ...
    'tpause', 0.01, ...
    'nTRskip', 10);

return;

