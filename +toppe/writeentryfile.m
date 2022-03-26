function writeentryfile(entryFile, varargin)
% function writeentryfile(entryFile, varargin)
%
% Convenience function for writing a toppe<CV1>.entry file.
% Ensures correct file format.

% defaults
arg.filePath       = '/usr/g/research/pulseq/myscan/';
arg.moduleListFile = 'modules.txt';
arg.loopFile       = 'scanloop.txt';
arg.b1ScalingFile  = 'tipdown.mod';
arg.readoutFile    = 'readout.mod';
arg.seqStampFile   = 'seqstamp.txt';

% substitute with provided keyword arguments
arg = toppe.utils.vararg_pair(arg, varargin);

% write to file
fid = fopen(entryFile, 'wt');
fprintf(fid, '%s\n', arg.filePath);
fprintf(fid, '%s\n', arg.moduleListFile);
fprintf(fid, '%s\n', arg.loopFile);
fprintf(fid, '%s\n', arg.b1ScalingFile);
fprintf(fid, '%s\n', arg.readoutFile);
fprintf(fid, '%s\n', arg.seqStampFile);
fclose(fid);


