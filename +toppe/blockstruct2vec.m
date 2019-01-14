function d = blockstruct2vec(block)
% Convert block struct to a [1 16] vector, suitable as a line entry in scanloop.txt (see loop2txt.m)
%
% Input:
%  block     struct containing TOPPE settings for the execution of one module
%
% Output:
%  d         [1 16] vector containing values to be written to scanloop.txt. See TOPPE manual and loop2txt.m for details.
%
% Example:  see examples/SPGR/main.m in https://github.com/toppeMRI/toppe

max_pg_iamp = 2^15-2;    % max amplitude in hardware units
dabon = 1;

% initialize
d = zeros(1,16);

% replace with those values that are provided
try
	d(1) = block.module;
catch
	error('module number must be specified');
end
if isfield(block, 'rfscale')
	d(2) = 2*round(max_pg_iamp * block.rfscale/2);   % force even amp
end
if isfield(block, 'thscale')
	d(3) = 2*round(max_pg_iamp * block.thscale/2);
else
	d(3) = max_pg_iamp;
end
if isfield(block, 'gxscale')
	d(4) = 2*round(max_pg_iamp * block.gxscale/2);
end
if isfield(block, 'gyscale')
	d(5) = 2*round(max_pg_iamp * block.gyscale/2);
end
if isfield(block, 'gzscale')
	d(6) = 2*round(max_pg_iamp * block.gzscale/2);
end
if isfield(block, 'dabslice')
	d(7) = block.dabslice;
end
if isfield(block, 'dabecho')
	d(8) = block.dabecho;
end
if isfield(block, 'dabview')
	d(9) = block.dabview;
end
if isfield(block, 'dabmode')
	d(10) = block.dabmode;
else
	d(10) = dabon;
end
if isfield(block, 'rot')
	d(11) = 2*round(max_pg_iamp * block.rot/2);
end
if isfield(block, 'rfphs')
	d(12) = 2*round(max_pg_iamp * block.rfphs/2);
end
if isfield(block, 'recphs')
	d(13) = 2*round(max_pg_iamp * block.recphs/2);
end
if isfield(block, 'textra')
	d(14) = block.textra;
end
if isfield(block, 'rffreq')
	d(15) = block.rffreq;
end
if isfield(block, 'wavnum')
	d(16) = block.wavnum;
else
	d(16) = 1;
end

return;


