function [gx,gy,gz,fname] = makeepi(fov, npix, zres, ncycles, varargin)
% function [gx,gy,gz] = makeepi(fov, npix, zres, ncycles, slthick, varargin)
%
% Not finished yet. Currently only does one echo (spin-warp).
%
% Make non-flyback multi-echo readout. Play out with TOPPE pulse sequence.
%
% INPUTS:
%  fov        - cm
%  npix       - number of pixels along fov
%  ncycles    - number of cycles of phase across slthick, added to x gradient 
%  zres       - resolution along z (partition) encoding direction (mm)
% Options:
%  mxg        Gauss/cm
%  mxs        Gauss/cm/ms
%  oprbw      receive bandwidth. Default: 125/4 kHz
%  slthick    cm
%
% $Id: makeepi.m,v 1.10 2018/11/13 18:41:39 jfnielse Exp $

import toppe.*
import toppe.utils.*
import toppe.utils.epi.*

% defaults
arg.mxg    = 4;    % Gauss/cm
arg.mxs    = 10;   % Gauss/cm/msec
arg.isdess = 0;
arg.system = systemspecs();
arg.flip   = false;
arg.oprbw  = 125/4;  % kHz
arg.separateMods = false;
arg.slthick = [];    % for writing spoiler to separate .mod file

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);      % requires MIRT

if arg.oprbw > 125
	error('oprbw can''t be larger than 125 kHz');
end

dt = 4e-3;             % gradient/daq sample duration (msec)

% set readout amplitude and number of datapoints to acquire, based on oprbw
gamma = 4.2575;                     % kHz/Gauss
decimation = 125/arg.oprbw;
g = (1/dt)/(gamma*fov)/decimation;             % Gauss/cm
npixro = npix*decimation;
if g > arg.mxg
	g = arg.mxg;
	fprintf(1, 'max g reduced to %.2f G/cm \n', g);
end

% readout gradient (ramp and plateu) 
gxro = g*ones(1,npixro);         
s = arg.mxs*dt;
ramp = 0:s:g ;  

% x prewinder.
area = sum([ramp gxro fliplr(ramp)])*dt*1e-3;                % G/cm*s
%gxprew = -trapwave(area/2, dt*1e-3, arg.mxg, arg.mxs*1e3);   % sqrt(3) since up to 3 gradients are playing simultaneously
gxprew = -trapwave2(area/2, arg.mxg, arg.mxs, dt);   % sqrt(3) since up to 3 gradients are playing simultaneously

% put together the whole waveform
gcrush = makecrusher(ncycles,fov/npix,0, arg.mxs, arg.mxg);
if ncycles > 0.5
	areacrush = sum(gcrush)*dt*1e-3;   % G/cm*sec
	arearo = sum([gxprew ramp gxro fliplr(ramp)])*dt*1e-3;
	bridge = mybridged(areacrush-arearo,gxro(end), arg.mxg, arg.mxs*1e3);  
	gx = [gxprew ramp gxro bridge fliplr(ramp)];
else
	gx = [ramp gxro fliplr(ramp)];
	area = sum(gx)*dt*1e-3;
	gbal = -trapwave2(area/2, arg.mxg, arg.mxs, dt);   
	gx = [gxprew ramp gxro fliplr(ramp) gbal];
end
if ncycles == 0.5
	gx = [gxprew ramp gxro fliplr(ramp)];
end

readstart = length([gxprew ramp]);
readend = readstart + length([gxro]);

iref = readstart + npixro/2;   % center of echo (point from which to calculate TE)

% create y and z phase-encodes
yarea = sum([gxro])*dt*1e-3;
gyprew = -trapwave2(yarea/2, arg.mxg, arg.mxs, dt);   
gy = 0*gx;
%gy((readstart-length(gyprew)-1):(readstart-2)) = gyprew;
gy(1:numel(gyprew)) = gyprew;
if ~arg.isdess
	gy((readend+1):(readend+length(gyprew))) = -gyprew;
end
if numel(gy)>numel(gx)
	gx = [gx zeros(1,numel(gy)-numel(gx))];
end

npixz = 30;
zfov = npixz*zres/10;   % cm
gzamp = (1/dt)/(gamma*zfov);             % Gauss/cm
gzro = gzamp*ones(1,npixz);
zarea = sum([gzro])*dt*1e-3;
gzprew = -trapwave2(zarea/2, arg.mxg, arg.mxs, dt);   
gz = 0*gx;
gz(1:numel(gzprew)) = gzprew;
if ~arg.isdess
	gz((readend+1):(readend+length(gzprew))) = -gzprew;
end
if numel(gz)>numel(gx)
	gx = [gx zeros(1,numel(gz)-numel(gx))];
	gy = [gy zeros(1,numel(gz)-numel(gy))];
end

%plot(gy); hold on; plot(gx,'r'); plot(gz,'g');

% let's make npre even
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
npad = mod(numel(gx),4);
gx = [gx; zeros(npad,1)];
gy = [gy; zeros(npad,1)];
gz = [gz; zeros(npad,1)];

if arg.flip
	gx = flipud(gx);
	gy = flipud(gy);
	gz = flipud(gz);
end

% write z phase-encode gradient to file (useful for spiral scans with toppe)
fname = sprintf('zphaseencode-zres%.1fmm-%s.mod', zres, date);
if arg.separateMods
	writemod('gz', gzprew(:), 'ofname', fname, 'desc', 'z phase-encode lobe');
end

% write to readout.mod
if arg.isdess==1
	fname = sprintf('gre-fov%dmm-npix%d-ncycles%.1f-oprbw%.3f-zres%dmm-isdess%d-flip%d-%s.mod', ...
	10*fov, npix, ncycles, arg.oprbw, zres, arg.isdess, arg.flip, date);
else
	fname = sprintf('gre-fov%dmm-npix%d-ncycles%.1f-oprbw%.3f-zres%dmm-%s.mod', ...
	10*fov, npix, ncycles, arg.oprbw, zres, date);
end
hdrints = [npre npixro iref];    % some useful numbers for recon
hdrfloats = [arg.oprbw];
writemod('gx', gx(:), 'gy', gy(:), 'gz', gz(:), 'ofname', fname, ...
         'desc', 'spin-warp (GRE) waveform', 'hdrints', hdrints, 'hdrfloats', hdrfloats);

% write spoiler
if arg.separateMods
	fname = sprintf('spoiler-zres%dmm-ncycles%d-%s.mod', arg.slthick*10,ncycles,date);
	gz = gcrush(:);
	gz = makecrusher(ncycles, arg.slthick, 0, arg.mxs/sqrt(2), arg.mxg/sqrt(2));
	writemod('gx', gz(:), 'gy', gz(:), 'gz', gz(:), 'ofname', fname, 'desc', 'spoiler');
end

return;

