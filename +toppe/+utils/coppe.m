% Packages toppe files into toppe-scanfiles.tgz, then copies it to the scanner
% Assumes you have SSH keys set up to log into romero
% Pretty much for UM fMRI lab use only...

fprintf('Making archive...');
[~,~] = system('tar czf toppe-scanfiles.tgz modules.txt scanloop.txt *.mod'); fprintf('done!\n');
fprintf('Copying to romero...');
[~,~] = system('scp toppe-scanfiles.tgz fmrilab@romero:~/amos/'); fprintf('done!\n');
fprintf('Copying to scanner...');
[~,~] = system('ssh -q fmrilab@romero "~/amos/pushtoppefiles"'); fprintf('done!\n');
disp('Ready to scan.');