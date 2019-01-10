3D SPGR sequence example

%% Set system limits.
% Good to back off a bit on slew to avoid PNS (even if scanner can slew at higher rate).
% Otherwise, accept default values.
sys = toppe.systemspecs('maxSlew', 15);  

%% make RF excitation .mod file

flip = 15;                  % flip angle (degrees)
slthick = 10;               % slab thickness (cm)
tbw = 8;                    % time-bandwidth product of SLR pulse 
dur = 2;                    % pulse duration (msec)
ncycles = 2*100;            % number of cycles of phase (spoiler) across slthick (0 = balanced)
ofname = 'tipdown.mod';     % output file name
ftype = 'min';              % minimum-phase SLR pulse (good for 3D imaging)
[rf, g] = toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncycles, 'ftype', ftype, 'ofname', ofname);

return;

writeloop;                     % create scanloop.txt

system('rm scan,spgr.tgz');
system('tar czf scan,spgr.tgz modules.txt scanloop.txt *.mod *.m ');

figure; toppe.plotseq(5000,5010); subplot(511); title('SPGR');
%playseq('scanloop.txt',3,0);


