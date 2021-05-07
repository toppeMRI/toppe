function coppe(varargin)
% function coppe(varargin)
%
% Load data for one echo (or all) from Pfile, EXCEPT dabslice=0 slot (which can contain corrupt data).
%
% Input options:
%  target          which scanner to copy pulse sequence to
%                  options: 'inside', 'outside'
%
% Packages toppe files into toppe-scanfiles.tgz, then copies it to the scanner
% Assumes you have SSH keys set up to log into romero/toro
% For UM fMRI lab use only...

import toppe.utils.*

arg.target        = 'inside'; % Default to inside scanner
arg = vararg_pair(arg, varargin);

fprintf('Making archive...');
[~,~] = system('tar czf toppe-scanfiles.tgz modules.txt scanloop.txt *.mod'); fprintf('done!\n');

try
    switch arg.target
        case 'inside'
            fprintf('Copying to romero...');
            [~,~] = system('scp toppe-scanfiles.tgz fmrilab@romero:~/amos/'); fprintf('done!\n');
            fprintf('Copying to scanner...');
            [~,~] = system('ssh -q fmrilab@romero "~/amos/pushtoppefiles"'); fprintf('done!\n');
            
        case 'outside'
            % switch to this version if we can debug firewall stuff (MWH note 4/12/21)
              fprintf('Copying to toro...');
              [~,~] = system('scp toppe-scanfiles.tgz fmrilab@toro:~/toppe_utils/'); fprintf('done!\n');
              fprintf('Copying to scanner...');
              [~,~] = system('ssh -q fmrilab@toro "~/toppe_utils/pushtoppefiles"'); fprintf('done!\n');
            
            % temp verion that requires the password to be entered
%             system('scp toppe-scanfiles.tgz fmrilab@toro:~/toppe_utils/');
%             system('ssh -q fmrilab@toro "~/toppe_utils/pushtoppefiles"');
        otherwise
            fprintf('Invalid target, valid targets are ''inside'' or ''outside''.\n');
    end
    disp('Files sucessfully copied. Finding # of slices...');
    scanloop_struc = importdata('scanloop.txt');
    scanloop_struc_data = scanloop_struc.data;
    max_sli = scanloop_struc_data(2);
    fprintf('Set # of slices on scanner to %d.\n',max_sli+1);
    disp('Ready to scan.');
catch
    fprintf('\n...failed. Not ready to scan.\n');
end