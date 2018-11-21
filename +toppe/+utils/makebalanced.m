function gbal = makebalanced(g, varargin)
% Add a trapezoid at end of g to make the total area zero.
%
% function gbal = makebalanced(g, varargin)
%
% Inputs:
%   g          G/cm
% Options:
%   maxSlew    G/cm/ms. Default: 10.
%   maxGrad    G/cm. Default: 5.
%   system     struct specifying hardware system limits, see systemspecs.m
%
% $Id: makebalanced.m,v 1.6 2018/10/26 01:38:45 jfnielse Exp $

import toppe.*
import toppe.utils.*

%% parse inputs
% Defaults
arg.maxSlew = 10;
arg.system  = toppe.systemspecs();
arg.maxGrad = arg.system.maxGrad;

%arg = toppe.utils.vararg_pair(arg, varargin);
arg = vararg_pair(arg, varargin);

dt = arg.system.raster;        % GE raster time (sec)

if ~exist('maxSlew', 'var')
	maxSlew = arg.system.maxSlew;
end
if ~exist('maxGrad', 'var')
	maxGrad = arg.system.maxGrad;
end

% ramp to zero
dg = -sign(g(end))*arg.maxSlew*dt*1e3;      % G/sample
g = [g(:)' g(end):dg:0];

area = sum(g)*dt;   % G/cm*sec
gblip = trapwave(abs(area), dt, arg.maxGrad, arg.maxSlew*1000);
gbal = [g(:); -sign(area)*gblip(:)];

return
