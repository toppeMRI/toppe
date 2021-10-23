function gcrush = makecrusher(ncycles,opslthick,system,gzarea,mxs,mxg)
% function gcrush = makecrusher(ncycles,opslthick,system,gzarea,mxs,mxg)
% 
% INPUTS:
%   ncycles     -- number of cycles of phase across slice/slab of thickness 'opslthick'
%   opslthick   -- cm
%   system      -- struct specifying hardware system info
%   gzarea      -- half-area of slice-select gradient [G/cm*sec]. Default: 0.
%   mxs         -- max slew rate [G/cm/msec]. Default: 10.
%   mxg         -- max gradient [G/cm]. Default: 4.
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
% University of Michigan 2019

import toppe.*
import toppe.utils.*

if ~exist('gzarea','var')
	gzarea = 0;
end
if ~exist('mxs','var')
	mxs = 10;   % G/cm/msec  % be very conservative about PNS
end
if ~exist('mxg','var')
	mxg = 4.0;   % G/cm
end

gamma = 4257.5;                              % Hz/Gauss
area = ncycles/(gamma*opslthick);           % G/cm*sec
%fprintf('makecrusher: ncycles = %d, opslthick = %.2f, area = %f\n', ncycles, opslthick, area);
%gcrush = trapwave(area-gzarea,4e-6,mxg,mxs*1e3);
dt = 4e-3;   % ms
gcrush = trapwave2(area-gzarea, mxg, mxs, dt);
gcrush = makeGElength(gcrush(:));

%% Check waveform against system hardware limits
if ~checkwaveforms(system, 'gz', gcrush)
	error('Waveforms failed system hardware checks -- exiting');
end

return;

% EOF
