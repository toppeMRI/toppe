function writeentryfile(entryFile, varargin)
% function writeentryfile(entryFile, varargin)
%
% Convenience function for writing a toppe<CV1>.entry file.
% Ensures correct file format.
%
% Optional keyword-argument inputs with defaults:
%  filePath       = '/usr/g/research/pulseq/myscan/';
%  version        = 7
%  moduleListFile = 'modules.txt';
%  loopFile       = 'scanloop.txt';
%  b1ScalingFile  = 'tipdown.mod';
%  readoutFile    = 'readout.mod';
%  seqStampFile   = 'seqstamp.txt';

% defaults
arg.filePath       = '/usr/g/research/pulseq/myscan/';
arg.version        = 7;
arg.moduleListFile = 'modules.txt';
arg.loopFile       = 'scanloop.txt';
arg.b1ScalingFile  = 'tipdown.mod';
arg.readoutFile    = 'readout.mod';
arg.seqStampFile   = 'seqstamp.txt';
arg.coresFile  = 'cores.txt';

% substitute with provided keyword arguments
arg = toppe.utils.vararg_pair(arg, varargin);

% write to file
fid = fopen(entryFile, 'wt');
if arg.version <= 6
    fprintf(fid, '%s\n', arg.filePath);
    fprintf(fid, '%s\n', arg.moduleListFile);
    fprintf(fid, '%s\n', arg.loopFile);
    fprintf(fid, '%s\n', arg.b1ScalingFile);
    fprintf(fid, '%s\n', arg.readoutFile);
    fprintf(fid, '%s\n', arg.seqStampFile);
    fprintf(fid, '%s\n', arg.coresFile);
else % version 7
    fprintf(fid, '%d\n', 1);
    fprintf(fid, '%s\n', arg.filePath);
end
fclose(fid);


