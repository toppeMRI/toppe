function [gx,gy,gz,fname] = makegre(fov, npix, zres, varargin)
%
% Make gradients for 3D spin-warp (cartesian gradient-echo) readout
%
% INPUTS:
%  fov        [1 1]    in-plane fov (cm)
%  npix       [1 1]    number of pixels along fov
%  zres       [1 1]    resolution along z (partition) encoding direction (cm)
% Options:
%  oprbw       Receive bandwidth. Default: +/- 125/4 kHz
%  ncycles     Number of cycles of phase across voxel dimension, added to x gradient. Default: 0 (balanced readout).
%  system      struct specifying system info, see systemspecs.m
%  ofname      Output file name. Default: 'readout.mod'
%  extrafiles  [bool] If true, writes z phase-encode and spoiler to separate .mod files.

zres = zres*10;   % mm

import toppe.*
import toppe.utils.*

% defaults
arg.oprbw  = 125/4;  % kHz
arg.ncycles = 0;
arg.system = systemspecs();
arg.ofname = 'readout.mod';
arg.extrafiles = false;
arg.flip   = false;
arg.isdess = 0;

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);      % requires MIRT

if arg.oprbw > 125
	error('oprbw can''t be larger than +/- 125 kHz');
end

% Gradient limits
mxg = arg.system.maxGrad;          
if strcmp(arg.system.gradUnit, 'mT/m')
	mxg = mxg/10;     % Gauss/cm
end
mxs = arg.system.maxSlew;
if strcmp(arg.system.slewUnit, 'T/m/s')
	mxs = mxs/10;     % Gauss/cm/ms
end

% gradient/daq sample duration (msec)
dt = 4e-3;             

% readout amplitude and number of datapoints to acquire
gamma = 4.2575;                     % kHz/Gauss
decimation = 125/arg.oprbw;
g = (1/dt)/(gamma*fov)/decimation;             % Gauss/cm
if g > mxg
	error(sprintf('Requested readout plateau gradient strength exceeds gradient amplitude limit (%.1f%%)', g/mxg*100));
end
npixro = npix*decimation;

% readout gradient
gxro = g*ones(1,npixro);        % plateau 
s = mxs*dt;                     % max change in gradient per sample
ramp = 0:s:g;

% x prewinder.
area = sum([ramp gxro fliplr(ramp)])*dt*1e-3;                % G/cm*s
%gxprew = -trapwave(area/2, dt*1e-3, mxg, mxs*1e3);   % sqrt(3) since up to 3 gradients are playing simultaneously
gxprew = -trapwave2(area/2, mxg, mxs, dt);   % sqrt(3) since up to 3 gradients are playing simultaneously

% put together the whole gx waveform
gcrush = makecrusher(arg.ncycles,fov/npix,0, mxs, mxg);
if arg.ncycles > 0
	areacrush = sum(gcrush)*dt*1e-3;   % G/cm*sec
	arearo = sum([gxprew ramp gxro fliplr(ramp)])*dt*1e-3;
	bridge = mybridged(areacrush-arearo,gxro(end), mxg, mxs*1e3);  
	gx = [gxprew ramp gxro bridge fliplr(ramp)];
else
	gx = [ramp gxro fliplr(ramp)];
	area = sum(gx)*dt*1e-3;
	gbal = -trapwave2(area/2, mxg, mxs, dt);   
	gx = [gxprew ramp gxro fliplr(ramp) gbal];
end

readstart = length([gxprew ramp]);
readend = readstart + length([gxro]);

iref = readstart + npixro/2;   % center of echo (point from which to calculate TE)

% create y phase encode
yarea = sum([gxro])*dt*1e-3;
gyprew = -trapwave2(yarea/2, mxg, mxs, dt);   
gy = 0*gx;
%gy((readstart-length(gyprew)-1):(readstart-2)) = gyprew;
gy(1:numel(gyprew)) = gyprew;
if ~arg.isdess
	gy((readend+1):(readend+length(gyprew))) = -gyprew;
end
if numel(gy)>numel(gx)
	gx = [gx zeros(1,numel(gy)-numel(gx))];
end

% create z (partition) encode
npixz = 30;             % only exists as helper variable used in next line
zfov = npixz*zres/10;   % cm
gzamp = (1/dt)/(gamma*zfov);             % Gauss/cm
gzro = gzamp*ones(1,npixz);
zarea = sum([gzro])*dt*1e-3;
gzprew = -trapwave2(zarea/2, mxg, mxs, dt);   
gz = 0*gx;
gz(1:numel(gzprew)) = gzprew;
if ~arg.isdess
	gz((readend+1):(readend+length(gzprew))) = -gzprew;
end
if numel(gz)>numel(gx)
	gx = [gx zeros(1,numel(gz)-numel(gx))];
	gy = [gy zeros(1,numel(gz)-numel(gy))];
end

% make npre even
npre = length([gxprew ramp]);
if mod(npre,2)
	gx = [0 gx];
	gy = [0 gy];
	gz = [0 gz];
	npre = npre+1;
	iref = iref+1;
end

gx = makeGElength(gx(:));
gy = makeGElength(gy(:));
gz = makeGElength(gz(:));

%[areacrush-sum(gx)*dt*1e-3]   % check that net gradient area is equal to areacrush 

% make number of waveform samples divisible by 4 (so that rhfrsize matches full length) 
%npad = mod(numel(gx),4);
%gx = [gx; zeros(npad,1)];
%gy = [gy; zeros(npad,1)];
%gz = [gz; zeros(npad,1)];

if arg.flip
	gx = flipud(gx);
	gy = flipud(gy);
	gz = flipud(gz);
end

%% Write to file(s)
%	fname = sprintf('gre-fov%dcm-npix%d-ncycles%.1f-oprbw%.3f-zres%dmm-%s.mod', ...
%	  fov, npix, arg.ncycles, arg.oprbw, 10*zres, date);
%end
hdrints = [npre npixro iref];    % some useful numbers for recon
hdrfloats = [arg.oprbw];
writemod('gx', gx(:), 'gy', gy(:), 'gz', gz(:), 'ofname', arg.ofname, ...
         'desc', '3D spin-warp (GRE) waveform', 'hdrints', hdrints, 'hdrfloats', hdrfloats);

if arg.extrafiles
	% write z phase-encode gradient to file (useful for spiral scans with toppe)
	%fname = sprintf('zencode-zres%.1fmm-%s.mod', zres, date);
	fname = 'zphaseencode.mod';
	writemod('gz', gzprew(:), 'ofname', fname, 'desc', 'z partition encoding lobe');

	% write z spoiler to separate .mod file
	if arg.ncycles > 0
		%fname = sprintf('spoiler-zres%dmm-ncycles%d-%s.mod', fov/npix,arg.ncycles,date);
		fname = 'spoiler.mod';
		gx = makecrusher(arg.ncycles, fov/npix, 0, mxs/sqrt(2), mxg/sqrt(2));
		gz = makecrusher(arg.ncycles, zres, 0, mxs/sqrt(2), mxg/sqrt(2));
		writemod('gx', gx(:), 'gz', gz(:), 'ofname', fname, 'desc', 'spoiler');
	end
end

return;

