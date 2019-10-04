function d = readloop(loopfile)
% read scanloop.txt (toppe.e scan loop definition, used with modules.txt).
%
% function d = readloop(loopfile)
%
% Input:
%    loopfile       default: 'scanloop.txt'

if ~exist('loopfile', 'var')
	loopfile = 'scanloop.txt';
end

% go past header
%fid = fopen(loopfile, 'r', 'ieee-be');
%s = fgets(fid);  % read a whole line
%s = fgets(fid);  % read a whole line
%nt       = fscanf(fid, '%d\t', 1);
%maxslice = fscanf(fid, '%d\t', 1);
%maxecho  = fscanf(fid, '%d\t', 1);
%maxview  = fscanf(fid, '%d\n', 1);
%fclose(fid);

% read data (tab-separated, starting at row 4)
d = dlmread(loopfile,'\t',3,0);

return;
