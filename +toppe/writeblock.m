function writeblock2file(ofname, blk, sys)
%
% Write a Pulseq block to a custome file format that
% can be read by the PulseGEq interpreter.

% open file for writing
fid = fopen(ofname, 'w', 'ieee-be');

% order is important
sub_writerf(fid, blk.rf, sys);
sub_writegrad(fid, blk.gx, sys);
sub_writegrad(fid, blk.gy, sys);
sub_writegrad(fid, blk.gz, sys);
sub_writeadc(fid, blk.adc, sys);

% done
fclose(fid);

return


%% Functions for writing channels to file
function sub_writerf(fid, rf, sys)

C = toppe.constants;

 % type
if isempty(rf)
    fwrite(fid, C.NULL, 'int16');
    return
end
fwrite(fid, C.ARBITRARY, 'int16');

% delay
fwrite(fid, round(rf.delay*1e6), 'int16');     % us

% amplitude
amp = max(abs(rf.signal/sys.gamma));     % Gauss
fprintf(fid, 'amp: %.5f\n', amp);

% waveform
rho = 2*round(abs(rf.signal/sys.gamma)/amp*C.max_pg_iamp/2);
theta = 2*round(angle(rf.signal)/pi*C.max_pg_iamp/2);
fwrite(fid, numel(rho), 'int16');   % number of samples in waveform
fwrite(fid, rho, 'int16');
fwrite(fid, theta, 'int16');

return


function sub_writegrad(fid, g, sys)

C = toppe.constants;

% type
if isempty(g)
    fwrite(fid, C.NULL, 'int16');
    return;
end
if strcmp(g.type, 'trap')
    fwrite(fid, C.TRAP, 'int16');
else
    fwrite(fid, C.ARBITRARY, 'int16');
end

% delay
fwrite(fid, round(g.delay*1e6), 'int16');     % us

% amplitude
amp = g.amplitude/sys.gamma/100;  % Gauss/cm
fprintf(fid, 'amp: %.5f\n', amp); 

% waveform
if strcmp(g.type, 'trap')
    fwrite(fid, round(g.riseTime*1e6), 'int16');
    fwrite(fid, round(g.flatTime*1e6), 'int16');
    fwrite(fid, round(g.fallTime*1e6), 'int16');
else
    % TODO
end

return


function sub_writeadc(fid, adc, sys)

C = toppe.constants;

if isempty(adc)
    fwrite(fid, C.NULL, 'int16');
    return
end
fwrite(fid, C.ADC, 'int16');

fwrite(fid, adc.numSamples, 'int16');   
fwrite(fid, round(adc.dwell*1e6), 'int16');   % us
fwrite(fid, round(adc.delay*1e6), 'int16');   % us

return







%% OLD CODE


%if ~checkwaveforms(system, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz)
%	('Waveforms failed system hardware checks -- exiting');
%end


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
%fprintf(fid, 'res:   %d\n', res);

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
