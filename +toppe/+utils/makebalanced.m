function gbal = makebalanced(g, varargin)
% Add a trapezoid at end of g to make the total area zero.
%
% function gbal = makebalanced(g, varargin)
%
% Inputs:
%   g          1D gradient waveform (G/cm)
% Options:
%   maxSlew    G/cm/ms. Default: 10.
%   maxGrad    G/cm. Default: 5.
%   system     struct specifying hardware system limits, see systemspecs.m
%              If 'system' is provided, it overrides 'maxSlew' and 'maxGrad' options.

import toppe.*
import toppe.utils.*

% parse inputs
arg.maxSlew = 10;
arg.maxGrad = 5;
arg.system  = []; 
arg = vararg_pair(arg, varargin);

if isempty(arg.system)
	maxSlew = arg.maxSlew;
	maxGrad = arg.maxGrad;
	sys = toppe.systemspecs;
	dt = sys.raster;        % GE raster time (sec)
else
	maxSlew = arg.system.maxSlew;  % assumes G/cm/ms
	maxGrad = arg.system.maxGrad;  % assumes G/cm
	dt = arg.system.raster;        % GE raster time (sec)
end

maxSlew = 0.995*maxSlew;    % so it passes hardware checks
maxGrad = 0.995*maxGrad;

% ramp to zero
dg = -sign(g(end))*maxSlew*dt*1e3;      % G/sample
g = [g(:)' g(end):dg:0];

% add balancing trapezoid
area = sum(g)*dt;    % G/cm*sec
gblip = trapwave2(abs(area), maxGrad, maxSlew, dt*1e3);
gbal = [g(:); -sign(area)*gblip(:)];

return
