function [gx, gy, gz] = makeepi(fov, N, nshots, sys, varargin)
% function [gx, gy (,gz)] = makeepi(fov, N, nshots, sys, varargin)
%
% Make single/multi-shot 2D EPI or 3D stack-of-EPI readout waveforms.
% Can be flyback or not; use ramp sampling or not.
%
% If writefiles = true, two .mod files are produced:
%   readout.mod     contains the echo train. Balanced gradients.
%   prephaser.mod   gx/gy(/gz) prephasing gradients (move to corner of kspace)
% An example of how these are used in a TOPPE sequence, 
% see TODO
%
% Required inputs:
%  fov        [1 2] or [1 3]   Field of view (cm)
%  N          [1 2] or [1 3]   Matrix size
%  nshots     number of RF shots needed to fully sample kx-ky plane
%  sys        struct specifying system info, see systemspecs.m
%
% Options:
%  Ry           (int) EPI undersampling factor. Default: 1
%  flyback      (true/false) Default: false
%  rampsamp     (true/false) Ramp sampling? Default: false
%  isbalanced   (true/false) Default: false
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
%
% Limitations/comments:
% This script does NOT check that echo spacings obey the scanner
% limits (to avoid mechanical resonances).

nx = N(1); ny = N(2);

if length(fov) == 3 & length(N) == 3
    is3d = true;
    nd = 3;
    nz = N(3);
    fovz = fov(3);
else
    is3d = false;
    nd = 2;
end

%% Defaults
arg.Ry = 1;
arg.flyback = false;
arg.rampsamp = false;
arg.isbalanced = false;
arg.decimation = 1;
arg.writefiles = false;

%% Substitute varargin values as appropriate and check inputs
arg = vararg_pair(arg, varargin);      % requires MIRT

if mod(nx,2)
	error('nx must be an even integer');
end

if mod(ny,2)
	error('ny must be an even integer');
end

if rem(arg.decimation, 1) | arg.decimation < 1
	error('decimation must be positive integer');
end

etl = ny/arg.Ry/nshots;

if mod(etl, 1)
    error('echo train length (= ny/Ry/nshots) must be an integer')
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

%% Some constants

dt = 4e-3;            % gradient/ADC sample duration (msec)
gamma = 4257.6;       % Hz/G

%% readout gradients

% y phase-encode blip
dky = arg.Ry*nshots/fov(2);    % ky spacing (cycles/cm)
gy.blip = toppe.utils.trapwave2(dky/(gamma), sys.maxGrad, sys.maxSlew/sqrt(2), dt);

% x readout gradient for one echo
if ~arg.rampsamp
    % no sampling on ramps
    gamp = 1/(fov(1)*gamma*dt*1e-3);    % Amplitude. Gauss/cm
    gx.plateau = gamp*ones(1,nx*arg.decimation);  % plateau 
    s = mxs*dt/sqrt(2);   % max change in gradient per sample
    gx.ramp = 0:s:gamp;

    % If needed, extend readout plateau to make room for y blips
    if length(gy.blip) > 2*length(gx.ramp)
        gx.ramp = linspace(0, gamp, ceil(length(gy.blip/2)));
    end

    % readout trapezoid (one echo in the EPI train)
    gx.echo = [gx.ramp gx.plateau fliplr(gx.ramp)];
else
    % sample on ramps
    % If non-flyback, make room for ky blips on turns
    res = fov(1)/nx;     % spatial resolution (cm)
    kmax = 1/(2*res);
    area = 2*kmax/gamma + (~arg.flyback)*dky/gamma;   % G/cm*s
    mxg = min(1/(fov(1)*gamma*dt*1e-3), sys.maxGrad);    % Must support Nyquist
    gx.echo = toppe.utils.trapwave2(area, mxg, sys.maxSlew/sqrt(2), dt);
end

% flyback gradient
if arg.flyback
    area = sum(gx.echo)*dt*1e-3;  % G/cm*s
    gx.flyback = -toppe.utils.trapwave2(area, sys.maxGrad, sys.maxSlew/sqrt(2), dt);
    if gx.flyback(end) == 0
        gx.flyback = gx.flyback(1:(end-1));
    end
end

% prewinders
area = sum(gx.echo)*dt*1e-3;  % G/cm*s
gx.pre = -toppe.utils.trapwave2(area/2, sys.maxGrad, sys.maxSlew/sqrt(nd), dt);

res = fov(2)/ny;     % cm
kmax = 1/(2*res);
area = 2*kmax/gamma;   % G/cm*s
gy.pre = -toppe.utils.trapwave2(area/2, sys.maxGrad, sys.maxSlew/sqrt(nd), dt);

if is3d
    res = fov(3)/nz;     % cm
    kmax = 1/(2*res);
    area = 2*kmax/gamma;   % G/cm*s
    gz.pre = -toppe.utils.trapwave2(area/2, sys.maxGrad, sys.maxSlew/sqrt(nd), dt);
else
    gz.pre = 0*gy.pre;
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

% convert to column vector and pad length to multiple of 4 samples
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

if arg.isbalanced
    % Make echo train balanced (zero gradient area).
    % This way the prephaser module has independent control over
    % the net gradient area.
    % Not the most efficient implementation though.

    % reduce slew a bit to reduce PNS
    area = sum(gy.et)*dt*1e-3;  % G/cm*s
    etrephase = -toppe.utils.trapwave2(area, sys.maxGrad, 0.7*sys.maxSlew/sqrt(2), dt);

    if arg.flyback
        iStop = length(gy.et) - length(gx.flyback); 
    else
        if ~arg.rampsamp
            iStop = length(gy.et) - length(gx.ramp);
        else
            iStop = length(gy.et) - floor(length(gy.blip)/2);
        end
    end

    gy.et = [gy.et(1:iStop); etrephase'];
    gy.et = [gy.et; zeros(length(gx.et)-length(gy.et),1)];
end

%% Assemble full waveform (mainly for viewing in Matlab)
gx.full = [gx.pre'; gx.et; -gx.pre'];
gy.full = [gy.pre'; gy.et; -gy.pre'];

if ~arg.writefiles
    return;
end


%% write to .mod files

gx.et = toppe.makeGElength(gx.et);
gy.et = toppe.makeGElength(gy.et);

% Store a few values in header for later use 
hdrints = [nx ny nshots arg.Ry length(gx.echo)]; 
if arg.flyback
    hdrints = [hdrints length(gx.flyback)];
end
if ~arg.rampsamp
    hdrints = [hdrints length(gx.ramp)];
end
hdrfloats = [fov(1) fov(2)];

toppe.writemod(sys, ...
    'gx', gx.et, 'gy', gy.et, ...
    'desc', 'EPI readout', ...
    'ofname', 'readout.mod', ...
    'hdrfloats', hdrfloats, ...
	'hdrints', hdrints);

toppe.writemod(sys, ...
    'gx', gx.pre', 'gy', gy.pre', 'gz', gz.pre', ...
    'desc', 'EPI prephasing gradients (move to corner of kspace)', ...
    'ofname', 'prephaser.mod');

return;



