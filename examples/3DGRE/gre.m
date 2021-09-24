function main
% 3D SPGR sequence example

% Set hardware limits (for design)
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
sys = toppe.systemspecs('maxSlew', 12.3, 'slewUnit', 'Gauss/cm/ms', ...
    'maxGrad', 5, 'gradUnit', 'Gauss/cm');


% Acquisition parameters
% fov and voxel size must be square (in-plane)
matrix = [240 240 10];
fov  = [24 24 10];       % cm
flip = 8;               % excitation flip angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

if matrix(1) ~= matrix(2) | fov(1) ~= fov(2)
    error('In-plane fov and voxel size must be square');
end

% Make .mod file containing z spoiler gradient and RF excitation
dur = 2;                    % RF pulse duration (msec)
slthick = fov(3)*0.8;       % Slab thickness (cm). A bit smaller than fov(3) to avoid aliasing.
ftype = 'min';              % minimum-phase SLR pulse (good for 3D imaging)
tbw = 8;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncyclesspoil*matrix(3), sys, ...
                       'ftype', ftype, 'ofname', 'tipdown.mod');

% Create readout.mod
toppe.utils.makegre(fov(1), matrix(1), fov(3)/matrix(3), sys, ... 
                    'ofname', 'readout.mod', 'ncycles', ncyclesspoil); 

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

% finalize scanloop.txt
toppe.write2loop('finish', sys);

% Play sequence in loop (movie) mode
%nModulesPerTR = 2;
%toppe.playseq(nModulesPerTR, 'nTRskip', 10);

% Create 'sequence stamp' file for TOPPE.
% This file is listed in the 5th row in toppe0.meta
toppe.preflightcheck('toppe0.meta', 'seqstamp.txt', sys);

% Write files to tar archive (for convenience only).
% toppe0.meta must be placed in /usr/g/research/pulseq/ on scanner host.
% The other files go in the folder specified on the first line in toppe0.meta.
system('tar cf gre.tar toppe0.meta modules.txt scanloop.txt *.mod seqstamp.txt');

return;

