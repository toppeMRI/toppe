function main
% 3D SPGR sequence example

%% Set system limits
% Good to back off a bit on slew to avoid PNS (even if scanner can slew at higher rate).
% Otherwise, accept default values.
sys = toppe.systemspecs('maxSlew', 10, 'maxGrad', 4);  

%% Acquisition parameters
% fov and voxel size must be square (in-plane)
matrix = [240 240 10];
fov  = [24 24 10];       % cm
flip = 10;               % excitation flip angle (degrees)
ncyclesspoil = 2;        % number of cycles of spoiler phase across voxel dimension (applied along x and z)

if matrix(1) ~= matrix(2) | fov(1) ~= fov(2)
	error('In-plane fov and voxel size must be square');
end

%% Make .mod file containing z spoiler gradient and RF excitation
ofname = 'tipdown.mod';     % Output file name
dur = 2;                    % RF pulse duration (msec)
slthick = fov(3)*0.8;       % Slab thickness (cm). A bit smaller than fov(3) to avoid aliasing.
ftype = 'min';              % minimum-phase SLR pulse (good for 3D imaging)
tbw = 8;                    % time-bandwidth product of SLR pulse 
toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncyclesspoil*matrix(3), ...
                       'ftype', ftype, 'ofname', ofname, 'system', sys);
%toppe.plotmod(ofname);

%% Create readout.mod
ofname = 'readout.mod';     % output file name
toppe.utils.makegre(fov(1), matrix(1), fov(3)/matrix(3), ... 
                    'system', sys, 'ofname', ofname, 'ncycles', ncyclesspoil); 
%toppe.plotmod(ofname);

%% Create scanloop.txt
rfmod = 1;           % module index, i.e., line number in modules.txt
readoutmod = 2;
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;    % degrees
ii = 1;                 % counts number of module executions
ny = matrix(2);
nz = matrix(3);

dabon = 1;
daboff = 0;
dabmode = dabon;   % best to leave acquisition on even during disdaqs so auto-prescan gets a signal (?)
waveform = 1;
textra = 0;        % add delay at end of module (int, microseconds)

loop = containers.Map;

for iz = 0:nz           % We'll use iz=0 for approach to steady-state
	for iy = 1:ny

		% rf excitation block (usage of 'block' here parallels its usage in Pulseq)
		%block = {'module', rfmod, ...
		%	'rfscale', 1.0, ...           % full amplitude (range is [-1 1]) 
		%	'gxscale', 0, ...
		%	'gyscale', 0, ...
		%	'gzscale', 1.0, ...
		%	'rfphs',   rfphs};
		block = [];
		block.module = rfmod;
		block.rfscale = 1.0;
		block.gxscale = 0;
		block.gyscale = 0;
		block.gzscale = 1.0;
		block.rfphs = angle(exp(1i*rfphs));   % radians
		d(ii,:) = sub_block2vec(block);
		ii = ii + 1;

		% readout
		block = [];
		block.module =  readoutmod;
		block.rfscale = 1.0;
		block.gxscale = 1.0;
		block.gyscale = ((iy-1+0.5)-ny/2)/(ny/2);    % phase-encode amplitude scaling, range is [-1 1]
		block.gzscale = ((iz-1+0.5)-nz/2)/(nz/2);    % partition-encode amplitude scaling
		block.dabslice = max(iz,0);                  % Convention: skip dabslice=0 
		block.dabecho = 0; 
		block.dabview = iy;                          % Convention: skip baseline (0) view
		block.recphs = angle(exp(1i*rfphs));         % radians
		d(ii,:) = sub_block2vec(block);
		ii = ii + 1;

		% update rf/rec phase
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt ;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
	end
end

% write loop array to scanloop.txt
toppe.loop2txt(d);

%% Play sequence in loop (movie) mode
%system('rm scan,spgr.tgz');
%system('tar czf scan,spgr.tgz modules.txt scanloop.txt *.mod *.m ');
%figure; toppe.plotseq(5000,5010); subplot(511); title('SPGR');
%playseq(2,0);

return;

% Convert paired cell array to a [1 16] line entry in scanloop.txt
function d = sub_block2vec(block)
% Input:
%  block     {'field 1', value1, 'field2', value2, ...} paired cell array containing TOPPE settings for the execution of one module
% Output:
% d          [1 16] vector containing values to be written to scanloop.txt. See TOPPE manual for details.
%            [tipdowncore ia_tipdown ia_th max_pg_iamp max_pg_iamp max_pg_iamp dabslice dabecho dabview daboff rot irfphase irfphase textra rffreq waveform]; 

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
	d(5) = 2*round(max_pg_iamp * block.gxscale/2);
end
if isfield(block, 'gzscale')
	d(6) = 2*round(max_pg_iamp * block.gxscale/2);
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


