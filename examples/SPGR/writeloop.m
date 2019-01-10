function d = writeloop
% 3D cartesian STFR scan loop, for TOPPE.
% $Id: writeloop.m,v 1.12 2019/01/09 20:17:15 jfnielse Exp $

load ../params   % ny nz freq_dz yresnom yres

[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = toppe.readmod('tipdown.mod');
nomflip = paramsfloat(11);  % flip angle of tipdown.mod at full amplitude (IV excitation pulse) (degrees)

% loop over tip angles
FLIPDOWN = [5:5:45] ;   % degrees
nscans = numel(FLIPDOWN);

dabon = 1;
daboff = 0;
max_pg_iamp = 2^15-2;  
ia_th = max_pg_iamp;

NL = 16;
nCoresPerTr = 3;
d = zeros(nCoresPerTr*ny*(nz+2)*(nscans),NL);

% loop through acquisitions
tipdowncore = 1;
daqcore     = 2;
spoilercore = 3;

textra = 0;
rffreq = freq_dz;

rf_spoil_seed = 117; %117;
rf_spoil_seed_cnt = 0;
rfphase = 0;              % radians

icnt = 1;
dabmode = dabon;   % best to leave acquisition on even during disdaqs so auto-prescan gets a signal (?)
waveform = 1;
textra = 0;

for iim = 1:nscans

	fprintf('%d ', iim);

	dabecho = iim-1;

	ia_tipdown = 2*round(FLIPDOWN(iim)/nomflip*max_pg_iamp/2);

	for iz = -1:nz   

		dabslice = max(iz,0);    % skip dabslice=0 as usual (disdaq)

		for iy = 1:ny
			if iz > 0
				ia_gy = 2*round(yresnom/yres* max_pg_iamp*(((iy-1+0.5)-ny/2)/(ny/2)) /2);   % y phase-encode amplitude
				dabview = iy;  % skip baseline (0) view as usual
				ia_gz = 2*round( max_pg_iamp*(((iz-1+0.5)-nz/2)/(nz/2)) /2);   % z phase-encode amplitude
			else   % disdaqs
				ia_gy = 0;
				dabview = 0;
				ia_gz = 0;
			end

			% calculate rf phase (rf spoiling)
			rfphase  = rfphase + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
			rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
			rfphasetmp = atan2(sin(rfphase), cos(rfphase));      % wrap phase to (-pi,pi) range
			irfphase = 2*round(rfphasetmp/pi*max_pg_iamp/2);     % even short int 

			textra = 0;

			d(icnt,:) = [tipdowncore ia_tipdown ia_th max_pg_iamp max_pg_iamp max_pg_iamp dabslice dabecho dabview daboff 0 irfphase irfphase textra rffreq waveform]; 
			icnt = icnt+1;

			d(icnt,:) = [daqcore 0 0 max_pg_iamp ia_gy ia_gz dabslice dabecho dabview dabmode 0 irfphase irfphase textra rffreq waveform];
			icnt = icnt+1;

			%textra = 20e3;
			d(icnt,:) = [spoilercore 0 0 max_pg_iamp max_pg_iamp max_pg_iamp dabslice dabecho dabview daboff 0 irfphase irfphase textra rffreq waveform];
			icnt = icnt+1;
		end
	end
end
fprintf('\n');

toppe.loop2txt(d(1:(icnt-1),:));        % write scanloop.txt

return;


