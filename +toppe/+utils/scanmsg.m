function scanmsg(entryFile)
% function scanmsg()
%
% Print scan instructions to terminal

fid = fopen(entryFile, 'r');
if fid == -1
    error('Failed to open entry file');
end
filePath = fgetl(fid);  % location of TOPPE scan files
fclose(fid);

fprintf('\nRename the .entry file and copy to /usr/g/research/pulseq/ on scanner host.\n');
fprintf('  (For example, name it ''toppe10.entry'' and set opuser1 (user CV1) = 10 when prescribing the scan.)\n');
fprintf(sprintf('Place all other files in %s on scanner host.\n', filePath));

