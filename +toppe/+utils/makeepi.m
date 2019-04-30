function [gx,gy] = makeepi(fov, N, nshots, varargin)
% function [gx,gy] = makeepi(fov, N, nshots, varargin)
%
% Make single/multi-shot EPI readout and write waveforms to a .mod file.
%
% INPUTS:
%  fov        cm (in-plane)
%  N          number of pixels (in-plane)
%  nshots     number of RF shots needed to fully sample (in-plane)
% Options:
%  ofname     default: 'readout.mod'
%  ncycles    number of cycles of phase across slthick, added at end of x and z gradient. Default: 0 (balanced readout).
%  system     struct specifying system info, see systemspecs.m

import toppe.*
import toppe.utils.*

%% Defaults
arg.ofname = 'readout.mod';
arg.oprbw  = 125;  % kHz
arg.ncycles = 2;
arg.system = systemspecs();
arg.zres = [];

%% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);      % requires MIRT

if mod(N,2)
	error('N must be an even integer');
end
if arg.oprbw > 125
	error('oprbw can''t be larger than (+/-) 125 kHz');
end
if rem(N, nshots)
	error('N/nshots must be an integer');
end

%% Gradient limits
mxg = 0.995*arg.system.maxGrad;          
if strcmp(arg.system.gradUnit, 'mT/m')
	mxg = mxg/10;     % Gauss/cm
end
mxs = 0.995*arg.system.maxSlew;
if strcmp(arg.system.slewUnit, 'T/m/s')
	mxs = mxs/10;     % Gauss/cm/ms
end

%% gradient/daq sample duration (msec)
dt = 4e-3;             

%% readout amplitude and number of datapoints to acquire
gamma = 4.2575;                     % kHz/Gauss
decimation = 125/arg.oprbw;
g = (1/dt)/(gamma*fov)/decimation;          % Gauss/cm
if g > mxg
	error(sprintf('Requested readout plateau gradient strength exceeds gradient amplitude limit (%.1f%%)', g/mxg*100));
end
npixro = N*decimation;    % number of 4us samples to acquire (per phase-encode)

%% readout gradient
plat = g*ones(1,npixro);        % plateau 
s = mxs*dt/sqrt(2);             % max change in gradient per sample. sqrt(2) since x and y gradients playing simultaneously
ramp = 0:s:g;
gro = [ramp plat fliplr(ramp)];
npre = length(ramp);            % number of samples before start of acquiring first echo (will be updated below)

%% x/y prewinders
area = sum(gro)*dt*1e-3;                % G/cm*s
gxpre = -trapwave2(area/2, mxg/sqrt(2), mxs/sqrt(2), dt);   % sqrt(n) since up to n gradients are playing simultaneously
area = sum(plat)*dt*1e-3;                % G/cm*s
gypre = -trapwave2(area/2, mxg/sqrt(2), mxs/sqrt(2), dt); 

%% EPI train
gx = [gro];                                                % first readout (echo)
gy = [0*gro];
ybliparea = area/N*nshots;
gyblip = trapwave2(ybliparea, mxg/sqrt(2), mxs/sqrt(2), dt);   % sqrt(2) since x and y gradients playing simultaneously
gyblip = [gyblip zeros(1,mod(length(gyblip),2))];               % make length even so we can divide by 2 next
nblip = length(gyblip);
etl = N/nshots;
for ii = 2:etl
	gx = [gx gro*(-1)^(ii-1)];
	gy = [gy(1:(end-nblip/2)) gyblip zeros(1,length(gro)-nblip/2)];
end

%% add pre-phaser gradient, and negate gx for every other shot
for ii = 1:nshots
	gxall(:,ii) = (-1)^(ii-1)*[gxpre(:); gx(:)];
	yscale = 1 - (2/N)*(ii-1);
	gyall(:,ii) = [zeros(length(gxpre)-length(gypre),1); gypre(:)*yscale; gy(:)];
end
npre = npre + length(gxpre);
	
%% add crusher to x 
l = 0;
for ii = 1:nshots
	gxarea = sum(gxall(:,ii))*dt*1e-3;                       % G/cm*s
	area = arg.ncycles/(1e3*gamma*fov/N) - gxarea;           % G/cm*s
	gtmp{ii} = trapwave2(area, mxg/sqrt(2), mxs/sqrt(2), dt)';
	%gcrush = makecrusher(arg.ncycles,fov/N,0, mxs/sqrt(2), mxg/sqrt(2));
	if length(gtmp{ii}) > l
		l = length(gtmp{ii});
	end
end
for ii = 1:nshots
	gcrush(:,ii) = [gtmp{ii}; zeros(l-length(gtmp{ii}),1)];
end
gxall = [gxall; gcrush]; %repmat(gcrush(:),1,size(gxall,2))];

%% add balancing gradient to y
%gyall = [gyall; repmat(0*gcrush(:),1,size(gyall,2))];   % make some space
l = 0;
for ii = 1:nshots
	area = sum(gyall(:,ii))*dt*1e-3;         % G/cm*s
	gtmp{ii} = -trapwave2(area, mxg/sqrt(2), mxs/sqrt(2), dt)';
	if length(gtmp{ii}) > l
		l = length(gtmp{ii});
	end
end
for ii = 1:nshots
	grew(:,ii) = [gtmp{ii}; zeros(l-length(gtmp{ii}),1)];
end
gyall = [gyall; grew];

%% write to .mod file
gx = makeGElength(gxall);
gy = makeGElength(gyall);
hdrints = [N nshots length(gro) npre]; 
hdrfloats = [fov];
writemod('ofname', arg.ofname, 'gx', gx, 'gy', gy, 'desc', 'EPI readout', ...
	'hdrints', hdrints);

return;



