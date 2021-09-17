function coppe(varargin)
% function coppe(varargin)
%  Function for **internal** University of Michigan fMRI lab use only to
%  send toppe files automatically to our scanners. Can be generalized for
%  other sites but the server names need to be updated.
%
%  A rough copy of the bash script called in this function (pushtoppefiles)
%  that we have on lab servers that can access the scanner is the 
%  following:
%
%       #!/bin/bash
%       # Support file for copying toppe scan files to the scanner, then unpacking them
%       set -e
%       IP=$(getscannerip)
%       echo "Copying files to scanner..."
%       ssh sdc@$IP 'scp fmrilab@toro:~/toppe_utils/toppe-scanfiles.tgz /usr/g/bin/; cd /usr/g/bin/; tar -xzf toppe-scanfiles.tgz'
%
%
% Input options:
%  target   which scanner to copy pulse sequence to, UM options: 'inside', 'outside'
%  use_pw   option to allow command line input, in the case a pw is needed
%
% Packages toppe files into toppe-scanfiles.tgz, then copies it to the scanner
% Assumes you have SSH keys set up to log into romero/toro

import toppe.utils.*

arg.target = 'inside'; % Default to inside scanner
arg.use_pw  = false;    % assume we have ssh keys setup to not need a pw
arg = vararg_pair(arg, varargin);

fprintf('Making archive...');
[status,cmdout] = system('tar czf toppe-scanfiles.tgz modules.txt scanloop.txt *.mod'); fprintf('done!\n');

if status
    error(cmdout)
end

try
    switch arg.target
        case 'inside'
            server_str = 'romero';
        case 'outside'
            server_str = 'toro';            
        otherwise
            fprintf('Invalid target, valid targets are ''inside'' or ''outside''.\n');
    end
    
    %% create linux commands with correct server target and run
    
    cmd1 = ['scp toppe-scanfiles.tgz fmrilab@',server_str,':~/toppe_utils/'];
    cmd2 = ['ssh -q fmrilab@',server_str,' "~/toppe_utils/pushtoppefiles"'];
    
    % 1. Send toppe files to server
    fprintf(['Copying to ',server_str,'...']);
    if ~arg.use_pw
        [st1, cmdo1] = system(cmd1);
    else
        [st1, cmdo1] = system(cmd1,'-echo');
    end
    if st1, error(cmdo1); end; fprintf('done!\n');
    
    % 2. Send toppe files from server to scanner computer and untar
    fprintf('Copying to scanner and untaring ...');
    if ~arg.use_pw
        [st2, cmdo2] = system(cmd2);
    else
        [st2, cmdo2] = system(cmd2,'-echo');
    end
    if st2, error(cmdo2); end; fprintf('done!\n');
    
    
    %% Find minimum number of slices
    disp('Files sucessfully copied. Finding # of slices...');
    scanloop_struc = importdata('scanloop.txt');
    scanloop_struc_data = scanloop_struc.data;
    max_sli = scanloop_struc_data(2);
    fprintf('Set minimum # of slices on scanner to %d.\n',max_sli+1);
    
    disp('Ready to scan.');
catch
    fprintf('\n...failed. Not ready to scan.\n');
end


