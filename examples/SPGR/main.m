% 3D SPGR sequence example

%import toppe.*
%import toppe.utils.*
%import toppe.utils.rf.*
%import toppe.utils.rf.jpauly.*

%% Set system limits
% Good to back off a bit on slew to avoid PNS (even if scanner can slew at higher rate).
% Otherwise, accept default values.
sys = toppe.systemspecs('maxSlew', 10, 'maxGrad', 4);  

%% Acquisition parameters
% fov and voxel size must be square (in-plane)
matrix = [240 240 50];
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
ii = 0;                 % counts number of module executions
ny = matrix(2);
nz = matrix(3);

dabmode = dabon;   % best to leave acquisition on even during disdaqs so auto-prescan gets a signal (?)
waveform = 1;
textra = 0;        % add delay at end of module (int, microseconds)

for iz = 0:nz           % We'll use iz=0 for approach to steady-state
	for iy = 1:ny

		% rf excitation
		ii = ii + 1;
		module = rfmod;   
		rfscale = 1.0;               % full amplitude (range is [-1 1]) 
		gxscale = 0;
		gyscale = 0;
		gzscale = 1.0;      %
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
		rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
		d(icnt,:) = [tipdowncore ia_tipdown ia_th max_pg_iamp max_pg_iamp max_pg_iamp dabslice dabecho dabview daboff 0 irfphase irfphase textra rffreq waveform]; 

		% readout
		ii = ii + 1;
		loop{ii}.module = readoutmod;   
		loop{ii}.gxscale = 1.0;
		loop{ii}.gyscale = ((iy-1+0.5)-ny/2)/(ny/2);   % phase-encode amplitude scaling
		if iz > 0
			loop{ii}.gzscale = ((iz-1+0.5)-nz/2)/(nz/2);   % partition-encode amplitude scaling
			loop{ii}.dabslice = max(iz,0);    % Convention: skip dabslice=0 
			loop{ii}.dabecho = 0;
			loop{ii}.dabview = iy;            % Convention: skip baseline (0) view
		else
			loop{ii}.gzscale = 0;    % approach to steady-state
			loop{ii}.dabslice = 0;
			loop{ii}.dabecho = 0;
			loop{ii}.dabview = 0;
		end
		loop{ii}.recphs = loop{ii}.rfphs;
	end
end

toppe.loopstruct2txt(loop);

return;

writeloop;                     % create scanloop.txt

system('rm scan,spgr.tgz');
system('tar czf scan,spgr.tgz modules.txt scanloop.txt *.mod *.m ');

figure; toppe.plotseq(5000,5010); subplot(511); title('SPGR');
%playseq('scanloop.txt',3,0);


