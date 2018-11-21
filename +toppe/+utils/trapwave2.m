function waveform = trapwave2(area, mxg, mxs, rasterTime) 
% Create gradient trapezoid with given area.
%
% Like trapwave.m, but avoids overshoots at plateau edges.
%
% function waveform = trapwave2(area, mxg, mxs, varargin) 
%
% Inputs:
%  area        G/cm*sec
%  mxg         G/cm
%  mxs         G/cm/ms
%  rasterTime  ms (On GE this should be 4e-3)
% 
% $Id: trapwave2.m,v 1.9 2018/11/15 14:26:55 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/trapwave2.m,v $

area = area*1e3;  % G/cm*msec

import toppe.*
import toppe.utils.*

% reduced peak amp/slew so it passes hardware checks
mxg = 0.99*mxg;
mxs = 0.99*mxs;

% construct trapezoid
dt = rasterTime;           % ms
tr = mxg/mxs;		         % Time for ramp to full scale.
Acrit = mxs*tr^2;	         % Area of maximum width triangle.
dg = mxs * dt;             % max change in gradient per raster time (G/cm)
if (area < Acrit) 
	rtime = sqrt(area/mxs);
	n = ceil(rtime/dt);
	ramp = 0:dg:(n*dg);
	waveform = [ramp fliplr(ramp)];
else
	nr = ceil(tr/dt);
	ramp = ([0:(nr-1)])*mxs*dt;
	areaRamps = 2*sum(ramp)*dt;
	np = ceil((area-areaRamps)/mxg/dt);
	plat = mxg*ones(1,np);
	waveform = [ramp plat hflip(ramp)];	
end

% scale down to desired area
wavArea = sum(waveform)*dt;           % Gauss/cm*ms
if wavArea < area
	error('Can''t scale down to desired area. Bug in code');
end
waveform = waveform/wavArea*area;

% waveforms must begin and end with zero (TOPPE convention)
waveform = [0  waveform 0];

% duration must be on 4 sample (16 us) boundary (TOPPE convention)
waveform = makeGElength(waveform(:))';

return;

