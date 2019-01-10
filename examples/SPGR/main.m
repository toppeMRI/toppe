% 3D SPGR sequence example

%import toppe.*
%import toppe.utils.*
%import toppe.utils.rf.*
%import toppe.utils.rf.jpauly.*

%% Set system limits
% Good to back off a bit on slew to avoid PNS (even if scanner can slew at higher rate).
% Otherwise, accept default values.
sys = toppe.systemspecs('maxSlew', 10, 'maxGrad', 4);  

%% Acquisition parameters
% fov and voxel size must be square (in-plane)
matrix = [240 240 50];
fov  = [24 24 10];       % cm
flip = 10;               % excitation flip angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

if fov(1) ~= fov(2)
	error('In-plane field-of-view must be square');
end

vox = fov./matrix;      % voxel dimensions (cm).

if vox(1) ~= vox(2)
	error('In-plane voxel size must be square');
end

%% Make .mod file containing z spoiler gradient and RF excitation
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
slthick = fov(3)*0.8;       % Slab thickness (cm). A bit smaller than fov(3) to avoid aliasing.
ftype = 'min';              % minimum-phase SLR pulse (good for 3D imaging)
tbw = 8;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncyclesspoil*matrix(3), ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
toppe.plotmod(ofname);

%% Create readout.mod
ofname = 'readout.mod';     % output file name
toppe.utils.makegre(fov(1), matrix(1), fov(3)/matrix(3), ... 
                    'system', sys, 'ofname', ofname, 'ncycles', ncyclesspoil); 
toppe.plotmod(ofname);

return;

writeloop;                     % create scanloop.txt

system('rm scan,spgr.tgz');
system('tar czf scan,spgr.tgz modules.txt scanloop.txt *.mod *.m ');

figure; toppe.plotseq(5000,5010); subplot(511); title('SPGR');
%playseq('scanloop.txt',3,0);


