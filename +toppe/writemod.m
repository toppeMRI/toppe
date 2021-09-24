function writemod(system, varargin)
% Write waveforms to .mod file, for use with toppe psd on GE scanners.
%
% function writemod(varargin)
%
% Assumes raster time (sample duration) of 4e-6 sec for all waveforms.
%
% Examples:
% >> writemod('rf', rho.*exp(1i*theta), 'gz', gzwaveform, 'nomflip', 30);
% >> writemod('gz', gzwaveform, 'desc', 'my spoiler gradient');
% >> lims = toppe.systemspecs('maxGrad', 130, 'gradUnit', 'mT/m');
% >> writemod('rf', myrf, 'gx', gzwav, 'system', lims);
%
% Input:
%   system        struct specifying hardware system info, see systemspecs.m
%
% Input options:
%   rf            Complex RF waveform, [ndat nrfpulses]
%   gx            [ndat ngxpulses]
%   gy            [ndat ngypulses]
%   gz            [ndat ngzpulses]
%   rfUnit        'Gauss' (default) or 'mT'
%   gradUnit      'Gauss/cm' (default) or 'mT/m'
%   ofname        Output filename.
%   desc          Text string (ASCII) descriptor.
%   nomflip       Excitation flip angle (degrees). Stored in .mod float header,
%                 but has no influence on b1 scaling. Default: 90.
%   hdrfloats     Additional floats to put in header (max 12)
%   hdrints       Additional ints to put in header (max 29)
%   nChop         [1 2] (int, multiple of 4) trim (chop) the start and end of
%                 the RF wavevorm (or ADC window) by this many 4us samples.
%                 Using non-zero nChop can reduce module duration on scanner.
%                 Default: [48 48] (= [192us 192us])

import toppe.*
import toppe.utils.*

nReservedInts = 2;   % [nChop(1) rfres], rfres = # samples in RF/ADC window
maxHdrInts = 32 - nReservedInts;

nChopDefault = [48 48];

%% parse inputs
% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];
arg.rfUnit    = 'Gauss';
arg.gradUnit  = 'Gauss/cm';
arg.ofname    = 'out.mod'; 
arg.desc      = 'TOPPE module';
arg.nomflip   = 90;
arg.hdrfloats = [];
arg.hdrints   = [];
arg.nChop = nChopDefault;

%arg = toppe.utils.vararg_pair(arg, varargin);
arg = vararg_pair(arg, varargin);

%% Check nChop
if arg.nChop(1) < nChopDefault(1) | arg.nChop(2) < nChopDefault(2)
    msg = ['nChop < 48 samples. Module duration (on scanner) will be ', ...
          'extended to account for RF/ADC dead/ringdown times as needed.'];
    warning(msg);
end

if ~isempty(arg.rf)
    if sum(abs(arg.rf([1:arg.nChop(1) (end-arg.nChop(2)+1):end]))) > 0
        msg = ['RF waveform must be zero during the first/last ', ...
              sprintf('%d/%d (nChop) samples.', arg.nChop(1), arg.nChop(2))]; 
        error(msg);
    end
end

if mod(arg.nChop(1),4) | mod(arg.nChop(2),4)
    error('Each entry in nChop must be multiple of 4');
end

%% Detect all-zero RF waveform, treat as empty, and give warning to user
if ~isempty(arg.rf) & norm(abs(arg.rf(:,1))) == 0
	arg.rf = [];
	warning('(First) RF waveform contains all zeros -- ignored');
end

%% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
	wavtype = fields{ii} ;   % 'rf', 'gx', 'gy', or 'gz'
	cmd = sprintf('%s = %s;', wavtype, sprintf('arg.%s', wavtype));
	eval(cmd);
end

%% Convert to Gauss and Gauss/cm
if strcmp(arg.rfUnit, 'mT')
	rf = rf/100;   % Gauss
end
if strcmp(arg.gradUnit, 'mT/m')
	gx = gx/10;    % Gauss/cm
	gy = gy/10;
	gz = gz/10;
end

%% Force all waveform arrays to have the same dimensions 
[rf, gx, gy, gz] = padwaveforms('rf', rf, 'gx', gx, 'gy', gy, 'gz', gz);
	
% Fixes to avoid idiosyncratic issues on scanner
%[rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp);

%% Check waveforms against system hardware limits
if ~checkwaveforms(system, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz)
	error('Waveforms failed system hardware checks -- exiting');
end

%% Header arrays
[paramsfloat] = sub_myrfstat(abs(rf(:,1,1)), arg.nomflip, system);
if ~isempty(arg.hdrfloats)
	if length(arg.hdrfloats) > 12
		error('max number of extra floats in .mod file header is 12');
	end
	paramsfloat(20:(19+length(arg.hdrfloats))) = arg.hdrfloats;  % populate header with custom floats 
end
%paramsint16 = [0 size(rf,1) arg.nChop]; % the first three ints are used by interpreter
paramsint16 = [arg.nChop(1) size(rf,1)-sum(arg.nChop)]; % the first three ints are used by interpreter
if ~isempty(arg.hdrints)
	if length(arg.hdrints) > maxHdrInts
		error(sprintf('max number of custom ints in .mod file header is %d', maxHdrInts));
	end
	paramsint16((nReservedInts+1):(nReservedInts+length(arg.hdrints))) = arg.hdrints;  % add custom ints
end

%% Write to .mod file
arg.desc = sprintf('Filename: %s\n%s', arg.ofname, arg.desc);
sub_writemod(arg.ofname, arg.desc, rf, gx, gy, gz, paramsint16, paramsfloat, system);

return;



function paramsfloat = sub_myrfstat(b1, nom_fa, system);
% Calculate RF parameters needed for RFPULSE struct in .e file.
% Needed for B1 scaling, SAR calculations, and enforcing duty cycle limits.
% See also mat2signa_krishna.m
%
% b1         real 1D vector containing B1 amplitude, size Nx1 [Gauss]
% nom_fa     nominal flip angle (degrees)

g = 1;  % legacy dummy value, ignore
nom_bw = 2000;

dt = 4e-6;                        % use 4 us RF sample width
%gamma = 4.2575e3;                  % Hz/Gauss
tbwdummy = 2;

hardpulse = max(abs(b1)) * ones(length(b1),1);    % hard pulse of equal duration and amplitude

pw            = length(b1)*dt*1e3;                                       % ms
nom_pw        = length(b1)*dt*1e6;                                        % us

if max(abs(b1)) == 0   % non-zero RF pulse
	error('RF waveform cannot be zero');
end
abswidth      = sum(abs(b1)) / sum(abs(hardpulse));
effwidth      = sum(b1.^2)   / sum(hardpulse.^2);
% or equivalently:  effwidth = sum(b1.^2)/(max(abs(b1))^2)/length(b1)
area          = abs(sum(b1)) / abs(sum(hardpulse)); 
dtycyc        = length(find(abs(b1)>0.2236*max(abs(b1)))) / length(b1);
maxpw         = dtycyc;
num           = 1;
max_b1        = system.maxRF;                       	% Gauss. Full instruction amplitude (32766) should produce max_b1 RF amplitude,
																		% as long as other RF .mod files (if any) use the same max_b1.
max_int_b1_sq = max( cumsum(abs(b1).^2)*dt*1e3 );   	% Gauss^2 - ms
max_rms_b1    = sqrt(mean(abs(b1).^2));              	% Gauss
nom_fa        = nom_fa;                              	% degrees
%nom_bw        = tbwdummy / (dt * length(b1));       	% Hz
% max_int_b1    = abs(sum(b1))*dt*1000

% calculate equivalent number of standard pulses
stdpw = 1;                                                % duration of standard pulse (ms)
stdpulse = 0.117 * ones(round(stdpw/(dt*1e3)),1);
numstdpulses = num * effwidth * (pw/stdpw) * (max(abs(b1))/0.117)^2;

%pulse50 = 0.117/180*50 * ones(round(stdpw/(dt*1e3)),1);
%num50 = sum(abs(b1).^2)/sum(pulse50.^2);

paramsfloat = [pw     abswidth  effwidth area          dtycyc      ...
              maxpw  num       max_b1   max_int_b1_sq max_rms_b1  ...
              90 nom_pw    nom_bw   g numstdpulses  nom_fa            ];  % hardcode 'opflip' to 90

%else % RF pulse is zero. This avoids division by zero.
%	paramsfloat = [pw 0 0 0 0 ...
%               0 1 0 0 ...
%               0 nom_pw    nom_bw   1 0];
%               0 
%end
               
return;


%%
function sub_writemod(fname,desc,rf,gx,gy,gz,paramsint16,paramsfloat,system)
%
%   desc           ASCII description 
%   rf             size(rho) = N x 1, where N = # waveform samples, 
%                  and Ncoils = number of transmit channels
%   gx,gy,gz       (Gauss/cm). vector of length N
%   paramsint16    Vector containing various int parameters, e.g., [npre ndur [additional optional values]]. nax length = 32.
%   paramsfloat    Vector containing various RF pulse parameters needed for SAR and B1 scaling
%                  calculations -- just a placeholder for now. Max length = 32.
%
% Jon-Fredrik Nielsen, jfnielse@umich.edu

rho = abs(rf);
theta = angle(rf);

npulses = size(rf,2);

% max length of params* vectors
nparamsint16 = 32;
nparamsfloat = 32;

% RF waveform is scaled relative to system.maxRF.
% This may be 0.25G/0.125G for quadrature/body RF coils (according to John Pauly RF class notes), but haven't verified...
if strcmp(system.rfUnit, 'mT')
	system.maxRF = system.maxRF/100;   % Gauss
end

% gradient waveforms are scaled relative to system.maxGrad
if strcmp(system.gradUnit, 'mT/m')
	system.maxGrad = system.maxGrad / 10;          % Gauss/cm
end

b1max = max(abs(rho(:)));     % Gauss

if (numel(paramsint16)>nparamsint16)
  error('writemod:nparamsint16',   'too many int16 parameters');
end
if (numel(paramsfloat)>nparamsfloat)
  error('writemod:nparamsfloat', 'too many float parameters');
end

% convert params* to row vectors, and pad to max length
paramsint16  = reshape(paramsint16,   1, numel(paramsint16));
paramsfloat  = reshape(paramsfloat, 1, numel(paramsfloat));
paramsint16  = [paramsint16   zeros(1, nparamsint16-numel(paramsint16))];
paramsfloat  = [paramsfloat zeros(1, nparamsfloat-numel(paramsfloat))];

[res,nrfpulses,ncoils] = size(rho);

fid = fopen(fname, 'w', 'ieee-be');

% peak gradient amplitude across all pulses
gmax = max(1.0, max(abs([gx(:); gy(:); gz(:)])));  % to avoid division by zero 

% write header
globaldesc = sprintf('RF waveform file for ssfpbanding project.\n');  
globaldesc = sprintf('%sCreated by %s.m on %s.\n', globaldesc, mfilename('fullpath'), datestr(now));  
fs = dbstack;
if numel(fs)>2
	globaldesc = sprintf('%scalled by (dbstack(3)): %s\n', globaldesc, fs(3).file);
end
globaldesc = sprintf('%sncoils = %d, res = %d\n', globaldesc, ncoils, res);  
desc = sprintf('%s%s\n', globaldesc, desc);
fwrite(fid, numel(desc), 'int16');      % number of characters in ASCII description
fwrite(fid, desc, 'uchar');

fwrite(fid, ncoils,  'int16');          % shorts must be written in binary -- otherwise it won't work on scanner 
fwrite(fid, res,     'int16');
fwrite(fid, npulses, 'int16');
fprintf(fid, 'b1max:  %f\n', system.maxRF);           % (floats are OK in ASCII on scanner)
fprintf(fid, 'gmax:   %f\n', gmax);

fwrite(fid, nparamsint16, 'int16');
fwrite(fid, paramsint16,  'int16');
fwrite(fid, nparamsfloat, 'int16');
for n = 1:nparamsfloat
	fprintf(fid, '%f\n', paramsfloat(n)); 
end

% write binary waveforms (*even* short integers -- the toppe driver/interpreter sets the EOS bit, so don't have to worry about it here)
max_pg_iamp = 2^15-2;                                   % RF amp is flipped if setting to 2^15 (as observed on scope), so subtract 2
rho   = 2*round(rho/system.maxRF*max_pg_iamp/2);
theta = 2*round(theta/pi*max_pg_iamp/2);
gx    = 2*round(gx/gmax*max_pg_iamp/2);
gy    = 2*round(gy/gmax*max_pg_iamp/2);
gz    = 2*round(gz/gmax*max_pg_iamp/2);

for ip = 1:npulses
	for ic = 1:ncoils
		fwrite(fid, rho(:,ip,ic), 'int16');
	end
	for ic = 1:ncoils
		fwrite(fid, theta(:,ip,ic),   'int16');
	end
	fwrite(fid, gx(:,ip),    'int16');
	fwrite(fid, gy(:,ip),    'int16');
	fwrite(fid, gz(:,ip),    'int16');
end

fclose(fid);

return;



%% function [rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp); %,tipupdelx,tipupdely)
function [rho,theta,gx,gy,gz] = sub_prepare_for_modfile(rho,theta,gx,gy,gz,addrframp); %,tipupdelx,tipupdely)

% make smooth RF ramp to avoid RF amp error
ramp = rho(1)*[linspace(0,1,5)]';
ramp = [-ramp/2; flipud(-ramp/2); ramp];
if ~addrframp
	ramp = [0; 0];   % to avoid non-zero gradients at beginning/end which causes problems with readmod
end

[n npulses ncoils] = size(rho);   % ncoils is number of RF transmit coils 
for ip = 1:npulses
	for ic = 1:ncoils
		rho2(:,ip,ic) = [ramp; rho(:,ip,ic)];
		theta2(:,ip,ic) = [0*ramp; theta(:,ip,ic)];
	end
end
rho = rho2;
theta = theta2;
gx = [zeros(numel(ramp),npulses); gx];
gy = [zeros(numel(ramp),npulses); gy];
gz = [zeros(numel(ramp),npulses); gz];


return;
