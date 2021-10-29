function [gx, gy] = makeepi(fov, N, nshots, sys, varargin)
% function [gx, gy] = makeepi(fov, N, nshots, sys, varargin)
%
% Make single/multi-shot 2D EPI readout and (optionally) 
% write waveforms to .mod file.
%
% Two .mod files are produced:
%   et.mod          contains the echo train. Balanced gradients.
%   prephaser.mod   gx/gy prephasing gradients (move to corner of kspace)
%
% Required inputs:
%  fov        [1 2] Field of view (cm)
%  N          [1 2] Matrix size
%  nshots     number of RF shots needed to fully sample
%  sys        struct specifying system info, see systemspecs.m
%
% Options:
%  Ry           (int) EPI undersampling factor. Default: 1
%  flyback      (true/false) Default: false
%  rampsamp     (true/false) Ramp sampling? Default: false
%  ncycles      (float) Number of cycles of phase across voxel width. Default: 2
%  decimation   (int) Design EPI readout as if ADC dwell time is 4us*decimation.
%               (The actual ADC dwell time is fixed to 4us in TOPPE)
%  writefiles   (true/false) Default: false
%
% Outputs:
%  gx        struct contain various elements of the waveform
%            gx.et = echo train portion
%            gx.pre = prephaser (written to prephaser.mod)
%            etc
%  gy        similar to gx
%
% Example:
%  >> sys = toppe.systemspecs();   % use default system values
%  >> [gx,gy] = toppe.utils.makeepi([24 20], [240 200], 40, sys, 'flyback', true);
%  >> plot([gx.et gy.et])

nx = N(1); ny = N(2);

%% Defaults
arg.Ry = 1;
arg.flyback = false;
arg.rampsamp = false;
arg.ncycles = 2;
arg.decimation = 1;
arg.writefiles = false;
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
mxg = 0.995*sys.maxGrad;          
if strcmp(sys.gradUnit, 'mT/m')
	mxg = mxg/10;     % Gauss/cm
end
mxs = 0.995*sys.maxSlew;
if strcmp(sys.slewUnit, 'T/m/s')
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
dky = arg.Ry*nshots/fov(2);    % ky spacing (cycles/cm)
gy.blip = toppe.utils.trapwave2(dky/(gamma), mxg, mxs/sqrt(2), dt);

% x readout gradient for one echo
if ~arg.rampsamp
    % no sampling on ramps
    g = (1/dt)/(gamma*1e-3*fov(1))/arg.decimation;  % amplitude. Gauss/cm
    gx.plateau = g*ones(1,nx*arg.decimation);  % plateau 
    s = mxs*dt/sqrt(2);   % max change in gradient per sample
    gx.ramp = 0:s:g;

    % If needed, extend readout plateau to make room for y blips
    if length(gy.blip) > 2*length(gx.ramp)
        gx.ramp = linspace(0, g, ceil(length(gy.blip/2)));
    end

    % readout trapezoid
    gx.echo = [gx.ramp gx.plateau fliplr(gx.ramp)];  % one echo in the EPI train
else
    % sample on ramps
    res = fov(1)/nx;     % spatial resolution (cm)
    kmax = 1/(2*res);    % cycles/cm
    gx.area = 2*kmax/gamma;   % G/cm * sec (area of each readout trapezoid)
    mxg = min(1/(fov(1)*gamma*dt*1e-3), sys.maxGrad)    % Gauss/cm
    gx.echo = toppe.utils.trapwave2(area, mxg, sys.maxSlew/sqrt(2), dt);

    % extend plateau to make room for ky blips on turns
    
end

% prewinders
gx.area = sum(gx.echo)*dt*1e-3;  % G/cm*s
gx.pre = -toppe.utils.trapwave2(gx.area/2, mxg, mxs/sqrt(3), dt);

if ~arg.rampsamp
    gy.area = sum(gx.plateau)*dt*1e-3;  % G/cm*s
else
    gy.area = sum(gx.echo)*dt*1e-3;  % G/cm*s
end
gy.pre = -toppe.utils.trapwave2(gy.area/2, mxg, mxs/sqrt(3), dt);

% flyback rewinder
if arg.flyback
    gx.flyback = -toppe.utils.trapwave2(gx.area, mxg, mxs/sqrt(2), dt);
    if gx.flyback(end) == 0
        gx.flyback = gx.flyback(1:(end-1));
    end
end

% assemble echo train
gx.et = [];
gy.et = [];
for iecho = 1:etl
    if arg.flyback
        gx.et = [gx.et gx.echo gx.flyback];
    else
        gx.et = [gx.et gx.echo*(-1)^(iecho+1)];
    end
end
gx.et = [gx.et 0];  % waveform must end with 0

% conver to column vector and pad length to multiple of 4 samples
gx.et = toppe.makeGElength(gx.et');

% add y blips
gy.et = 0*gx.et;
nb = length(gy.blip);
for ii = 1:(etl-1)
    if arg.flyback
        n = length(gx.echo) + length(gx.flyback);
        n2 = length(gx.echo) + floor(length(gx.flyback)/2);
        iStart = (ii-1)*n + n2 - nb/2;
    else
        n = length(gx.echo);
        iStart = ii*n - nb/2;
    end
    iStop = iStart + nb - 1;
    gy.et(iStart:iStop) = gy.blip;
end

% Add y rephaser (of etl portion only)
% so the etl portion is balanced.
% This way the prephaser module has independent control over
% the net gradient area.
% Not the most efficient implementation though.
gy.area = sum(gy.et)*dt*1e-3;  % G/cm*s
gy.etrephase = -toppe.utils.trapwave2(gy.area, mxg, mxs/sqrt(3), dt);
iStop = length(gy.et) - length(gx.ramp);
gy.et = [gy.et(1:iStop); gy.etrephase'];
gx.et = [gx.et; zeros(length(gy.et)-length(gx.et),1)];

if ~arg.writefiles
    return;
end

%% write to .mod files
% Store a few values in header that are useful for recon
gx.et = toppe.makeGElength(gx.et);
gy.et = toppe.makeGElength(gy.et);

hdrints = [nx ny nshots arg.Ry length(gx.ramp) length(gx.echo)]; 
if arg.flyback
    hdrints = [hdrints length(gx.flyback)];
end
hdrfloats = [fov(1) fov(2)];
toppe.writemod(sys, ...
    'gx', gx.et, 'gy', gy.et, ...
    'desc', 'EPI readout', ...
    'ofname', 'et.mod', ...
    'hdrfloats', hdrfloats, ...
	'hdrints', hdrints);

toppe.writemod(sys, ...
    'gx', gx.pre', 'gy', gy.pre', ...
    'ofname', 'prephaser.mod');

return;



