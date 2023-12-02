function [rf, g, freq, fnamestem] = makeslr(flip, slthick, tbw, dur, ncycles, system, varargin)
% Create slice-selective SLR pulse with gradient crusher (or balancing blip) before it.
%
% This code has gotten pretty messy and needs a cleanup for readability.
%
% function [rf, g, freq, fnamestem] = makeslr(flip, slthick, tbw, dur, ncycles, system, varargin)
% 
% Inputs:
%   flip                 flip angle (degrees)
%   slthick              slab thickness (cm)
%   tbw                  time-bandwidth product of SLR pulse (e.g., 4 or higher)
%   dur                  pulse duration (msec)
%   ncycles              number of cycles of phase (spoiler) (0 = balanced)
%   system               struct specifying hardware system info, see systemspecs.m
% Options:
%   ofname               output file name
%   type                 'st' (default), 'se', 'ex', ... (John Pauly's SLR toolbox: dzrf.m)
%   ftype                'ls' (default), 'min', ... 
%   sliceOffset          (cm) Default: 0. Determines the return value 'freq', to be used in scanloop.txt.
%   forBlochSiegert      Writes RF excitation and rephaser gradient into separate .mod files
%                        (rephaser goes after Fermi pulse)
%   writeModFile         true (default) or false
%   discardPrephaser     Don't include the balanced (pre-phaser) gradient trapezoid (if ncycles=0)? Default: false.
%   spoilDerate          [1 1], range is [0.1 1.0]. Derate slew rate by this factor during pre/rewinders. Default: 1.0.
%   doDisplay            boolean
% Outputs:
%   rf               Gauss
%   g                Gauss/cm
%   freq             RF frequency offset (for slice selection)
%   fnamestem
%
% Example:
%  >> sys = systemspecs('maxGrad', 3, 'maxSlew', 10);
%  >> [rf,g] = toppe.utils.rf.makeslr(180, 0.5, 8, 8, 4, 'type', 'se', 'system', sys);
%

import toppe.*
import toppe.utils.*
import toppe.utils.rf.*
import toppe.utils.rf.jpauly.*

if strcmp(flip, 'test')
	sub_test();
	return;
end

%% parse inputs
% Default values 
arg.ofname          = [];
arg.sliceOffset     = 0;
arg.type            = 'st';
arg.ftype           = 'ls';
arg.forBlochSiegert = false;
arg.writeModFile    = true;
arg.isPresto        = false;
arg.spoilDerate     = 1.0;
arg.doDisplay       = false;
arg.maxSlewScale    = 1;  % scales the max slew (useful for PNS control)
if arg.spoilDerate < 0.1 | arg.spoilDerate > 1.0
	error('spoilDerate must be in the range [0.1 1.0]');
end
if strcmp(arg.type, 'se')
	arg.discardPrephaser = true;
else
	arg.discardPrephaser = false;
end

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

if round(tbw) ~= tbw
%	error('tbw must be an integer');
end

mxg = system.maxGrad;          
mxs = system.maxSlew;
mxs = arg.maxSlewScale * mxs;

sliceOffset = arg.sliceOffset;

if ncycles==0 | arg.isPresto
	isBalanced=1;
else
	isBalanced=0;
end

dt = system.raster*1e-3;        % msec

% 
%switch arg.type
%	case {'sat'}
%		flip = 90;
%	case {'se','inv'}
%		flip = 180;
%end

%% Design rf pulse and slice-select gradient
resex = round(dur/dt);  % number of (4us) samples in RF waveform
[rfex,gex,irep,iref,gplateau] = sub_myslrrf(flip, dt*resex, tbw, arg.type, slthick, mxg, 0.99*mxs, arg.ftype, isBalanced, system, arg.spoilDerate);

% remove balancing (pre-phaser) gradient at beginning of pulse
if arg.discardPrephaser & ncycles == 0
	I = find(gex==0);
	J = find(I>10);      % go past any zeros that may be present at start of gradient
	gex = gex(I(J(1)):end);
	rfex = rfex(I(J(1)):end);
	gex = makeGElength(gex);
	rfex  = makeGElength(rfex);
end

% add spoiler gradient along slice-select
if ncycles > 0 & ~strcmp(arg.type, 'se')
	gspoil = toppe.utils.makecrusher(ncycles, slthick, system, 0, 0.99*mxs, mxg); 
	areaSpoil = sum(gspoil)*system.raster;  % G/cm*sec
	areass = sum(gex)*system.raster;        % area of slice-select gradient [G/cm*sec]
	if areass > areaSpoil
  		% no need to add spoiler (already big enough)
	else
		I = find(gex ~= 0);
		gdiff = diff(gex(I(1):end));
		J = find(gdiff==0);
		plateauStart = I(1)+J(1);
		ramp = [0; gex(1:(plateauStart-1))];
		gex = gex((plateauStart-1):end);
		rfex = rfex((plateauStart-1):end);
		bridge = toppe.utils.mybridged(areaSpoil-areass, gex(1), mxg, arg.spoilDerate*0.99*mxs*1e3); 
		npre = length([ramp; bridge']);
		gex = [ramp; bridge'; gex];
		rfex = [0*ramp; 0*bridge'; rfex];
	end
end

% make sure waveforms start and end at zero
rfex = [0; rfex; 0];
gex = [0; gex; 0];

% output file name
switch arg.type
	case 'ex'
		%rfex = rfex/90*flip;
		fnamestem = sprintf('tipdown-flip%d-slthick%.1fcm-tbw%.1f-dur%.1fms-ncycles%.1f-%s', flip, slthick, tbw, dur, ncycles, date);
	case 'se'
		fnamestem = sprintf('spinecho-slthick%.1fcm-tbw%.1f-dur%.1fms-ncycles%.1f-%s', slthick, tbw, dur, ncycles, date);
end

if ~isempty(arg.ofname)
	I = strfind(arg.ofname, '.mod');
	if ~isempty(I)
		fnamestem = arg.ofname(1:(I-1));
	else
		fnamestem = arg.ofname;   % '.mod' will be appended to output file name below
	end
end

% slice offset frequency
freq = system.gamma*gplateau*sliceOffset; % Hz

% crusher
if isBalanced | arg.isPresto
	gcrush = [];
else
	gcrush = makecrusher(ncycles,slthick,system,0,arg.spoilDerate*0.99*mxs,mxg);
end
gcrush = [gcrush(:); zeros(4,1)];

gspoil = [gcrush(:); 0*gex(:)];

%gex = [gcrush(:); gex(:) ];
%rfex = [0*gcrush(:); rfex(:) ];
%iref = iref + length(gcrush);
%irep = irep + length(gcrush);
switch arg.type
	case 'se'
		gex = [gcrush(:); gex; gcrush(:)];
		rfex = [0*gcrush(:);rfex; 0*gcrush(:)];
		gspoil = [gspoil; gcrush(:)];
end

if arg.isPresto
	gpresto1 = makecrusher(ncycles,slthick,system,0, 0.99*mxs, mxg);     % played before readout
	gpresto2 = makecrusher(2*ncycles,slthick,system,0, 0.99*mxs, mxg);   % played at end of TR
	gex = [gpresto2(:); gex(:); -gpresto1(:)];
	rfex = [0*gpresto2(:); rfex(:); 0*gpresto1(:)];
end

if arg.writeModFile
	rfex = toppe.makeGElength(rfex);
	gex = toppe.makeGElength(gex);
	%writemod('rf', rfex(:), 'gx', gspoil(:), 'gy', gspoil(:), 'gz', gex(:), ...
	writemod(system, 'rf', rfex(:), 'gz', gex(:), ...
		'nomflip', flip, 'ofname', sprintf('%s.mod',fnamestem), 'desc', 'SLR pulse');
end
%plotmod(sprintf('%s.mod',fnamestem));

% write RF trapezoid and rephaser to separate .mod files (for Bloch-Siegert with toppe.e)
if arg.forBlochSiegert
	rfex1 = [0; rfex(1:irep); 0];
	gex1  = [0; gex(1:irep);  0];
	rfex1 = makeGElength(rfex1);
	gex1 = makeGElength(gex1);
	if arg.writeModFile
		writemod(system, 'rf', rfex1(:), 'gz', gex1(:), 'nomflip', flip, ...
			'ofname', sprintf('%s-norep.mod',fnamestem), 'desc', 'SLR pulse, less rephaser');
	end
	rfex2 = [0; rfex((irep+1):end); 0];
	rfex2 = makeGElength(rfex2);
	gex2 = [0; gex((irep+1):end); 0];
	gex2 = makeGElength(gex2);
	if arg.writeModFile
		writemod(system, 'gz', gex2(:), 'nomflip', flip, ...
			'ofname', sprintf('%s-rephaser.mod',fnamestem), 'desc', 'SLR pulse, rephaser');
	end
end

% output
% make equal length
rf = rfex;
g  = gex; 
n = max(length(rf),length(g));
rf = [rf; zeros(n-length(rf),1)];
g = [g; zeros(n-length(g),1)];

% display slice profile
if arg.doDisplay
	Z = linspace(-slthick,slthick,100);
	m = slicesim([0 0 1],rf,g,dt,Z,1000,100,true);
end

return;



%% create RF pulse with slice-select gradient
function [rf,gex,irep,iref,gplateau,areaprep,idep,arearep] = sub_myslrrf(flip, dur, tbw, type, slthick, mxg, mxs, ftype, isBalanced, sys, spoilDerate)
%
% INPUTS:
%  flip         flip angle (degrees)
%  dur        - RF pulse duration (msec)
%  tbw        - total # zero crossings = time-bandwidth product
%  type       - 'st'/'ex'/'se' for small-tip / pi/2 excitation / spin-echo pulse
%  slthick    - cm 
%  mxg        G/cm
%  mxs        G/cm/ms
%  ftype      'ls' (default)
%  isBalanced default: true
%  sys        system struct  (systemspecs.m)
%
% OUTPUTS:
%  rf         - b1 waveform (real) (Gauss)
%  gss        - slice-select gradient waveform (Gauss/cm)
%  irep       - beginning of rephaser gradient (sample #)
%  iref       - location of "center" of RF pulse (point from which TE is to be calculated) (sample #)
%  gplateau   - amplitude of slice-select gradient (needed for off-isocenter slices)
%
% $Id: makeslr.m,v 1.27 2018/11/15 14:26:50 jfnielse Exp $

import toppe.*
import toppe.utils.*
import toppe.utils.rf.*
import toppe.utils.rf.jpauly.*

%% make rf waveform
dt = sys.raster*1e-3;                 % sample size (msec)
res = round(dur/dt);
dur = res*dt;

nrf = 200;                               % number of samples for SLR design
if round(tbw) ~= tbw & strcmp(type, 'st')
	error('tbw must be an integer for small-tip design');
end
rf = dzrf(nrf,tbw,type,ftype);       % row vector
rf = real(rf);
rf = resample(rf, res, nrf);

% scale to Gauss
rf = flip/180*pi * rf / sum(rf);     % normalized such that flip angle (radians) = sum(rf)
gamma = sys.gamma*1e-7;              % kHz/Gauss
rf = rf / gamma / dt /2/pi;          % Gauss

% make duration even
rf = [rf(:); zeros(mod(length(rf),2),1)];
npix = length(rf);

% find center (peak) of RF pulse
I = find(abs(rf) > max(abs(rf(:)))-eps);
iref= round(mean(I));

%% make slice-select gradient waveform 
bw = tbw / dur;                    % kHz
gplateau = bw / (gamma * slthick);     % slice-select gradient amplitude (G/cm)

if gplateau > mxg
	error('gradient exceeds mxg');
end

% slice-select trapezoid 
gss = gplateau*ones(1,npix);   % plateau of slice-select gradient
s = mxs * dt * 0.995;   % max change in g per sample (G/cm), slightly decreased to avoid floating point error
gss_ramp = [s:s:gplateau];
if isempty(gss_ramp)
	gss_ramp = 0;
end

if gplateau-gss_ramp(end) > s % Fix the boundary of ramp and plateau when g is not a multiple of s
    gss_ramp = [gss_ramp (gplateau+gss_ramp(end))/2];
end

gss_trap = [gss_ramp gss fliplr(gss_ramp)]; 
iref = iref + numel(gss_ramp);

% slice-select rephaser gradient
switch type
	case {'ex', 'st', 'sat'}
		arearep = (gss_trap(iref)/2+sum(gss_trap((iref+1):end))) * dt * 1e-3;            % G/cm*s
		gzrep = -trapwave2(arearep, mxg, spoilDerate*mxs, dt);
	case 'se'
		gzrep = [];
	case 'inv'
		gzrep = [];
end

% slice-select prephaser trapezoid
switch type
	case {'ex', 'st', 'sat'}
		%area = sum(gss( ((length(ramp)+1):(length(ramp)+midpoint)):end)) * dt * 1e-3             % G/cm*s
		areaprep = sum([gss_trap gzrep]) * dt * 1e-3;            % G/cm*s
		gzprep = -trapwave2(areaprep, mxg, spoilDerate*mxs, dt);
	case 'se'
		gzprep = [];
	case 'inv'
		gzprep = [];
end

% put together the full slice-select gradient and rf waveform
% make gss and rf the same length
if ~isBalanced
	gzprep = [];
end

irep = length([gzprep gss_trap])+1;
iref = iref + numel(gzprep);
gex = [gzprep gss_trap gzrep];
idep = numel(gzprep);

rf = [0*gzprep(:); zeros(length(gss_ramp),1); rf; zeros(length(gss_ramp)+length(gzrep),1)];

% rescale gradient rephaser to ensure accurate refocusing (flat phase across slice)
if strcmp(type, 'ex')
	Z = linspace(-0.5*slthick/2,0.5*slthick/2,50);
	m = slicesim([0 0 1],rf,gex,dt,Z,1000,100,false);
	ph = angle(m);  % we want this to be flat
	P = polyfit(Z,ph,1);
	gamma = 4.2576e3;   % Hz/Gauss
	extraarea = P(1)/(2*pi*gamma);    % G/cm*s
	gzrep = -trapwave2(arearep-extraarea, mxg, spoilDerate*mxs, dt);
	gex = [gzprep gss_trap gzrep];
	%m = slicesim([0 0 1],rf,gex,dt,Z,1000,100,true);
end

% ensure that duration is on a 16 us (4 sample) boundary
rf = makeGElength(rf(:));
gex = makeGElength(gex(:));

return;


function sub_test

% design a spin-echo pulse and simulate SE slice profile
flip = 180;
slthick = 2;   % cm
tbw = 6;
dur = 6;       % ms
ncycles = 8;   % number of cycles of spoiling across slthick
ofname = 'se.mod';
sys = toppe.systemspecs('maxSlew', 15, 'maxRf', 0.25);
[rf,gex,freq180] = toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncycles,...
	'system', sys, 'sliceOffset', 0, 'ofname', ofname, 'type', 'se');
toppe.plotmod(ofname);

m0 = [1 0 0];  % initial magnetization
dt = 4e-3;     % sample (raster, dwell) time (us)
Z = linspace(-slthick, slthick);
T1 = 1000; T2 = 100;
figure;
toppe.utils.rf.slicesim(m0,rf,gex,dt,Z,T1,T2);

% design a small-tip pulse and simulate slice profile
flip = 20;
slthick = 0.5;   % cm
tbw = 6;
dur = 2;       % ms
ncycles = 2;   % number of cycles of spoiling across slthick
ofname = 'st.mod';
[rf,gex] = toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncycles,...
	'system', sys, 'sliceOffset', 0, 'ofname', ofname, 'type', 'st');
toppe.plotmod(ofname);

m0 = [0 0 1];  % initial magnetization
Z = linspace(-slthick, slthick);
figure;
toppe.utils.rf.slicesim(m0,rf,gex,dt,Z,T1,T2);

return;
