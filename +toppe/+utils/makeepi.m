function [gx,gy] = makeepi(fov, N, nshots, varargin)
% function [gx,gy] = makeepi(fov, N, nshots, varargin)
%
% Make single/multi-shot 2D EPI readout and (optionally) write waveforms to a .mod file.
%
% Required inputs:
%  fov        [1 2] Field of view (cm)
%  N          [1 2] Matrix size
%  nshots     number of RF shots needed to fully sample
%
% Options:
%  Ry           (int) EPI undersampling factor. Default: 1
%  flyback      (true/false) Default: false
%  rampsamp     (true/false) Ramp sampling? Default: false
%  ncycles      (float) Number of cycles of phase across voxel width, added at end of x and z gradient. Default: 2
%  decimation   (int) Design EPI readout as if ADC dwell time is 4us*decimation.
%               (The actual ADC dwell time is still fixed to 4us in TOPPE. May change in future.)
%  system       struct specifying system info, see systemspecs.m
%  writemod     (true/false) Default: false
%  ofname       Output file name. Default: 'readout.mod' 

import toppe.*
import toppe.utils.*

nx = N(1); ny = N(2);

%% Defaults
arg.Ry = 1;
arg.flyback = false;
arg.rampsamp = false;
arg.ncycles = 2;
arg.decimation = 1;
arg.system = systemspecs();
arg.writemod = false;
arg.ofname = 'readout.mod';

%% Substitute varargin values as appropriate and check inputs
arg = vararg_pair(arg, varargin);      % requires MIRT

if mod(nx,2)
	error('nx must be an even integer');
end
if mod(ny,2)
	error('ny must be an even integer');
end
if rem(ny, nshots)
	error('ny/nshots must be an integer');
end
if rem(arg.decimation, 1) | arg.decimation < 1
	error('decimation must be positive integer');
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

%% gradient/ADC sample duration (msec)
dt = 4e-3;             

%% kspace steps
gamma = 4257.6;        % Hz/G

res = fov(1)/nx;       % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;   % G/cm * sec (area of each readout trapezoid)
dkx = 1/fov(1);

%dkz = 1/fov(3);        % kz spacing (cycles/cm) corresponding to Delta=1

etl = ny/arg.Ry/nshots;

if mod(etl, 1)
    error('echo train length (= ny/Ry/nshots) must be an integer')
end

%% readout gradients

% y phase-encode blip
dky = arg.Ry/fov(2);    % ky spacing (cycles/cm)
gy.blip = toppe.utils.trapwave2(dky/(gamma*1e3), mxg, mxs/sqrt(2), dt);

% x readout gradient for one echo
if ~arg.rampsamp
    gamma = 4.2575;              % kHz/Gauss
    g = (1/dt)/(gamma*fov(1))/arg.decimation;  % amplitdue. Gauss/cm
    gx.plateau = g*ones(1,nx*arg.decimation);  % plateau 
    s = mxs*dt/sqrt(2);   % max change in gradient per sample
    gx.ramp = 0:s:g;

    % If needed, extend readout plateau to make room for y blips
    if length(gy.blip) > 2*length(gx.ramp)
        gx.ramp = linspace(0, g, ceil(length(gy.bglip/2)));
    end

    % readout trapezoid
    gx.echo = [gx.ramp gx.plateau fliplr(gx.ramp)];  % one echo in the EPI train
else
    % ramp sampling
end

% x prewinder
gx.area = sum(gx.echo)*dt*1e-3;  % G/cm*s
gx.pre = -trapwave2(gx.area/2, mxg, mxs/sqrt(3), dt);

keyboard

%% x/y prewinders
area = sum(gro)*dt*1e-3;                % G/cm*s
gxpre = -trapwave2(area/2, mxg/sqrt(2), mxs/sqrt(2), dt); 
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
	%gcrush = makecrusher(arg.ncycles,fov/N,arg.system,0, mxs/sqrt(2), mxg/sqrt(2));
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
writemod(arg.system, ...
    'gx', gx, 'gy', gy, ...
    'desc', 'EPI readout', ...
    'ofname', arg.ofname, ...
	'hdrints', hdrints);

return;



