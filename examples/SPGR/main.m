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
rfmodnum = 1;           % module index, i.e., line number in modules.txt
readoutmodnum = 2;
rfphs = 0;              % radians
rf_spoil_seed_cnt = 0;
rf_spoil_seed = 117;      % degrees
for iz = 0:matrix(3)    % We'll use iz=0 for approach to steady-state
	for iy = 1:matrix(2)
		% rf excitation
		loop{ii}.module = rfmod;   
		loop{ii}.flip = 10;        % can't exceed design flip angle ('flip')
		loop{ii}.gxscale = 0;
		loop{ii}.gyscale = 0;
		loop{ii}.gzscale = 1;      % full amplitude (range is [-1 1])
		rfphs = rfphs + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
		loop{ii}.rfphs == angle(exp(1i*(rfphs)));

      ia_gz = 2*round( max_pg_iamp*(((iz-1+0.5)-nz/2)/(nz/2)) /2);   % z phase-encode amplitude
      else   % discarded acquisitions (disdaqs)
         ia_gz = 0;
      end
      dabslice = max(iz, 0);

      for il = 1:nl
         dabview = il;    % skip baseline (0) view as usual
         phi = round(2*pi*(il-1)/nl/pi*max_pg_iamp);

         % phase for rf spoiling
         rfphase = rfphase + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
         rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
         rfphasetmp = atan2(sin(rfphase), cos(rfphase));      % wrap phase to (-pi,pi) range
         irfphase = 2*round(rfphasetmp/pi*max_pg_iamp/2);     % short int 

	end
end



scanloop

return;

writeloop;                     % create scanloop.txt

system('rm scan,spgr.tgz');
system('tar czf scan,spgr.tgz modules.txt scanloop.txt *.mod *.m ');

figure; toppe.plotseq(5000,5010); subplot(511); title('SPGR');
%playseq('scanloop.txt',3,0);


