function [rf, g, freq, fnamestem] = makeslr(flip, slthick, tbw, dur, ncycles, varargin)
% Create slice-selective SLR pulse with gradient crusher (or balancing blip) before it.
%
% function [rf, g, freq, fnamestem] = makeslr(flip, slthick, tbw, dur, ncycles, varargin)
% 
% Inputs:
%   flip                 flip angle (degrees)
%   slthick              slab thickness (cm)
%   tbw                  time-bandwidth product of SLR pulse (e.g., 4 or higher)
%   dur                  pulse duration (msec)
%   ncycles              number of cycles of phase (spoiler) (0 = balanced)
% Options:
%   ofname               output file name
%   system               struct specifying hardware system info, see systemspecs.m
%   type                 'ex' (default), 'se', ... (John Pauly's SLR toolbox: dzrf.m)
%   ftype                'ls' (default), 'min', ... 
%   sliceOffset          (cm) Default: 0. Determines the return value 'freq', to be used in scanloop.txt.
%   forBlochSiegert      Writes RF excitation and rephaser gradient into separate .mod files
%                        (rephaser goes after Fermi pulse)
%   writeModFile         true (default) or false
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
% $Id: makeslr.m,v 1.27 2018/11/15 14:26:50 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+rf/makeslr.m,v $

import toppe.*
import toppe.utils.*
import toppe.utils.rf.*
import toppe.utils.rf.jpauly.*

%% parse inputs
% Default values 
arg.ofname          = [];
arg.sliceOffset     = 0;
arg.system          = toppe.systemspecs();
arg.type            = 'ex';
arg.ftype           = 'ls';
arg.forBlochSiegert = 0;
arg.writeModFile    = true;
arg.isPresto        = false;
arg.dorfmask        = false;

% Substitute varargin values as appropriate
arg = vararg_pair(arg, varargin);

if round(tbw) ~= tbw
	error('tbw must be an integer');
end


%% Design pulse and gradients

% multiply by 0.99 so gradients will pass checkwaveforms()
mxg = 0.99*arg.system.maxGrad;          
if strcmp(arg.system.gradUnit, 'mT/m')
	mxg = mxg/10;     % Gauss/cm
end
mxs = 0.99*arg.system.maxSlew;
if strcmp(arg.system.slewUnit, 'T/m/s')
	mxs = mxs/10;     % Gauss/cm
end

dorfmask    = arg.dorfmask;
sliceOffset = arg.sliceOffset;

if ncycles==0 | arg.isPresto
	isBalanced=1;
else
	isBalanced=0;
end

%tbw = 4;
dt = arg.system.raster*1e3;        % msec

% Design rf pulse and slice-select gradient
resex = round(dur/dt);  % number of (4us) samples in RF waveform
[rfex,gex,irep,iref] = sub_myslrrf(dt*resex, tbw, arg.type, slthick, mxg, mxs, arg.ftype, isBalanced, arg.system);
                               
% output file name
switch arg.type
	case 'ex'
		rfex = rfex/90*flip;
		fnamestem = sprintf('tipdown-flip%d-slthick%.1fcm-tbw%d-dur%.1fms-ncycles%.1f', flip, slthick, tbw, dur, ncycles);
	case 'se'
		fnamestem = sprintf('spinecho-slthick%.1fcm-dur%.1fms-ncycles%.1f-tbw%d-%s', slthick, dur, ncycles, tbw, date);
end

if ~isempty(arg.ofname)
	I = strfind(arg.ofname, '.mod');
	if ~isempty(I)
		fnamestem = arg.ofname(1:(I-1));
	else
		fnamestem = arg.ofname;   % '.mod' will be appended to output file name below
	end
end

% Bloch-Siegert modules
if arg.forBlochSiegert
	if arg.writeModFile
		writemod('rf', rfex(:), 'gz', gex(:), 'nomflip', flip, 'ofname', ...
		sprintf('%s-plateauOnly.mod',fnamestem), 'desc', 'SLR pulse', 'system', arg.system);
	end
end

% slice offset frequency
gamma = 4.257e3;  % Hz/Gauss
gplateau = gex(end/2);
freq = arg.system.gamma*gplateau*sliceOffset; % Hz

% crusher
if isBalanced | arg.isPresto
	gcrush = [];
else
	gcrush = makecrusher(ncycles,slthick,0,mxs,mxg);
end
gcrush = [gcrush(:); zeros(4,1)];

gspoil = [gcrush(:); 0*gex(:)];

gex = [gcrush(:); gex(:) ];
rfex = [0*gcrush(:); rfex(:) ];
iref = iref + length(gcrush);
irep = irep + length(gcrush);
switch arg.type
	case 'se'
		gex = [gex; gcrush(:)];
		rfex = [rfex; 0*gcrush(:)];
		gspoil = [gspoil; gcrush(:)];
end

if arg.isPresto
	gpresto1 = makecrusher(ncycles,slthick,0, mxs, mxg);     % played before readout
	gpresto2 = makecrusher(2*ncycles,slthick,0, mxs, mxg);   % played at end of TR
	gex = [gpresto2(:); gex(:); -gpresto1(:)];
	rfex = [0*gpresto2(:); rfex(:); 0*gpresto1(:)];
end

if arg.writeModFile
	writemod('rf', rfex(:), 'gx', gspoil(:), 'gy', gspoil(:), 'gz', gex(:), ...
		'nomflip', flip, 'ofname', sprintf('%s.mod',fnamestem), 'desc', 'SLR pulse', 'system', arg.system);
end
%plotmod(sprintf('%s.mod',fnamestem));

% write RF trapezoid and rephaser to separate .mod files (for Bloch-Siegert with toppe.e)
if arg.forBlochSiegert
	rfex1 = rfex(1:irep);
	gex1 = gex(1:irep);
	if arg.writeModFile
		writemod('rf', rfex1(:), 'gz', gex1(:), 'nomflip', flip, ...
			'ofname', sprintf('%s-norep.mod',fnamestem), 'desc', 'SLR pulse, less rephaser', 'system', arg.system);
	end
	rfex2 = rfex((irep+1):end);
	gex2 = gex((irep+1):end);
	if arg.writeModFile
		writemod('gz', gex2(:), 'nomflip', flip, ...
			'ofname', sprintf('%s-rephaser.mod',fnamestem), 'desc', 'SLR pulse, rephaser', 'system', arg.system);
	end
end

% output
rf = rfex;
g  = gex; 

return;



%% create RF pulse with slice-select gradient
function [rf,gss,irep,iref,gplateau,areaprep,idep,arearep] = sub_myslrrf(dur, tbw, type, slthick, mxg, mxs, ftype, isBalanced, sys)
% function [rf,gss,irep,iref,gplateau] = myslrrf(dur,tbw,type,slthick,isBalanced,type,ftype)
%
% INPUTS:
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
dt = sys.raster*1e3;                 % sample size (msec)
res = round(dur/dt);
dur = res*dt;

nrf = 200;                               % number of samples for SLR design
rf = dzrf(nrf,tbw,type,ftype);       % row vector
rf = real(rf);
rf = resample(rf, res, nrf);

% scale to Gauss
switch type
	case {'ex','sat'};
		flip = pi/2;
	case {'se','inv'}
		flip = pi;
end
rf = flip * rf / sum(rf);               % normalized such that flip angle (radians) = sum(rf)
gamma = sys.gamma*1e-3;               % kHz/Gauss
rf = rf / gamma / dt /2/pi;             % Gauss

% make duration even
rf = [rf(:); zeros(mod(length(rf),2),1)];
npix = length(rf);

iref = find(rf==max(rf));              % center of RF pulse

%% make slice-select gradient waveform 
bw = tbw / dur;                    % kHz
g = bw / (gamma * slthick);      % slice-select gradient amplitude (G/cm)
gplateau = g;

if g > mxg
	error('gradient exceeds mxg');
end

% slice-select trapezoid 
gss = g*ones(1,npix);                             % plateau of slice-select gradient
s = mxs * dt;                                     % max change in g per sample (G/cm)
ramp = [0:s:g];
gss = [ramp gss fliplr(ramp)]; 
iref = iref + numel(ramp);

% slice-select rephaser gradient
switch type
	case {'ex', 'st', 'sat'}
		%midpoint = find(rf==max(rf(:)));                          % center of main lobe
		%area = sum(gss((midpoint+1):end)) * dt * 1e-3  ;          % G/cm*s
		arearep = sum([gss((iref+1):end)]) * dt * 1e-3;            % G/cm*s
		gzrep = -trapwave2(arearep, mxg, mxs, dt);
		%gzrep = [gzrep zeros(1,mod(length(gzrep),2)) ];            % make length even
	case 'se'
		gzrep = [];
	case 'inv'
		gzrep = [];
end

% slice-select prephaser trapezoid
switch type
	case {'ex', 'st', 'sat'}
		%area = sum(gss( ((length(ramp)+1):(length(ramp)+midpoint)):end)) * dt * 1e-3             % G/cm*s
		areaprep = sum([gss(1:iref)]) * dt * 1e-3;            % G/cm*s
		gzprep = -trapwave2(areaprep, mxg, mxs, dt);
		gzprep = [gzprep zeros(1,mod(length(gzprep),2)) ];            % make length even
	case 'se'
		gzprep = [];
	case 'inv'
		gzprep = [];
end

% put together the full slice-select gradient
if ~isBalanced
	gzprep = [];
end

irep = length([gzprep gss]);
iref = iref + numel(gzprep);
gss = [gzprep gss gzrep];
idep = numel(gzprep);

% make gss and rf the same length. 
rf = [0*gzprep(:); zeros(length(ramp),1); rf; zeros(length(ramp)+length(gzrep),1)];

% ensure that duration is on a 16 us (4 sample) boundary
rf = makeGElength(rf(:));
gss = makeGElength(gss(:));

return;
