function coppe(varargin)
% function coppe(varargin)
%   Function for **internal** University of Michigan fMRI lab use only to
%   send toppe files automatically to our scanners. Can be generalized for
%   other sites but the server names need to be updated.
%
% Input options:
%   cv      CV number, used to determing where on the scanner to put files
%   target  which scanner to copy pulse sequence to, UM options: 'inside', 'outside'
%   force   option to force overwrite of existing sequences corresponding
%           to cv (default is 0, no overwrite)
%   user    user name on server to use for file transfer (defaults to current user name)
%   use_pw  option to allow command line input, in the case a pw is needed
%           (default is 0)
%   dir     target directory to write to under /srv/nfs/psd/usr/psd on scanner
%           (default is [username]/toppe[cv #])
%   version toppe software version (5 for tv5, 6 for tv6)
%
% Packages toppe files into toppe-scanfiles.tgz, then copies it to the scanner
% Assumes you have the following SSH keys set up:
%   1. ssh key from host --> server
%   2. ssh key from scanner --> host
%
% if SSH keys require a passcode, set use_pw = 1 to show prompt
%

    import toppe.utils.*
    basedir = '/srv/nfs/psd/usr/psd'; % base directory on scanner
    
    % Set default arguments
    arg.target = 'inside';
    arg.force = 0;
    [~,arg.user] = system('echo -n $USER');
    arg.use_pw  = 0;
    arg.cv = [];
    arg.dir = 'auto';
    arg.version = 6; %Toppe version
    arg = vararg_pair(arg, varargin);

    % set the default target directory
    if strcmpi(arg.dir,'auto')
        arg.dir = sprintf('%s/toppe%d', arg.user, arg.cv);
    end

    % get the host IP address
    [~,host_IP] = system('wget -qO- ifconfig.me/ip');
    
    % zip the sequence using tar
    fprintf('Zipping toppe files using tar...');
    if(arg.version < 6)
        [status,cmdout] = system('tar czf toppe-scanfiles.tgz modules.txt seqstamp.txt scanloop.txt *.mod');
    else
        [status,cmdout] = system('tar czf toppe-scanfiles.tgz modules.txt seqstamp.txt scanloop.txt cores.txt toppeN.entry *.mod');
    end
    if status
        error(cmdout)
    else
        fprintf(' SUCCESS\n');
    end

    % set server names
    switch arg.target
        case 'inside'
            server_str = 'epyc';
        case 'outside'
            server_str = 'goliath';
        otherwise
            error('Invalid target, valid targets are ''inside'' or ''outside''.\n');
    end
    
    % Set up the commands to run on the scanner
    cmd_scanner = [];
    if ~arg.force
        % check for existing entry file
        cmd_scanner = sprintf('%s if [[ -f %s/pulseq/toppe%d.entry ]]; then exit 15; fi;', ...
            cmd_scanner, basedir, arg.cv);
    end
    % check if directory exists, if it doesn't, create it
    cmd_scanner = sprintf('%s if [ ! -d %s/%s ]; then mkdir -p %s/%s; fi;', ...
            cmd_scanner, basedir, arg.dir, basedir, arg.dir);
    % cd into the directory
    cmd_scanner = sprintf('%s if ! cd %s/%s; then exit 16; fi;', ...
        cmd_scanner, basedir, arg.dir);
    % scp toppe files from host
    cmd_scanner = sprintf('%s if ! scp -q %s@%s:%s/toppe-scanfiles.tgz ./; then exit 17; fi;', ...
        cmd_scanner, arg.user, host_IP, pwd);
    % unzip the tar file
    cmd_scanner = sprintf('%s if ! tar -xzf toppe-scanfiles.tgz; then exit 18; fi; rm toppe-scanfiles.tgz;', ...
        cmd_scanner);
    % rename the entry file and move it to pulseq directory
    cmd_scanner = sprintf('%s if ! mv toppeN.entry %s/pulseq/v%d/toppe%d.entry; then exit 19; fi;', ...
        cmd_scanner, basedir, arg.version,arg.cv);
    % replace the first line with the right directory
    cmd_scanner = sprintf('%s if ! sed -i "1s#.*#%s/%s/#" %s/pulseq/v%d/toppe%d.entry; then exit 20; fi;', ...
        cmd_scanner, basedir, arg.dir, basedir, arg.version, arg.cv);

    % ssh into server, then into scanner, and run the scanner commands using bash
    cmd_host = sprintf('ssh -q %s@%s ssh -q sdc@10.0.1.1 /bin/bash << EOF\n%s\nEOF', ...
        arg.user, server_str, cmd_scanner);
    
    % use -echo to show the bash output (in case pw is required)
    if arg.use_pw
        eval_args = {cmd_host, '-echo'};
    else
        eval_args = {cmd_host};
    end

    % Copy the file to the scanner and unzip
    fprintf('Copying to scanner and unzipping...');
    [out,msg] = system(eval_args{:});
    switch out
        case 15
            error('cv number %d is already taken, use force argument to force overwrite', arg.cv);
        case 16
            error('failed to cd into %s/%s\noutput:\n%s', basedir, arg.dir, msg)
        case 17
            error('failed to scp from host to scanner, are your keys set up right?\noutput:\n%s', msg)
        case 18
            error('failed to unzip the tar file\noutput:\n%s', msg)
        case 19
            error('failed to move the entry file\noutput:\n%s', msg)
        case 0
            fprintf(' SUCCESS\nsequence files written to %s/%s\n', basedir, arg.dir);
        otherwise
            error('unknown error\noutput:\n%s',msg)
    end

end


