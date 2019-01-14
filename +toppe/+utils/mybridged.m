function g = mybridged(area_target,g0,mxg,mxs)
% function g = mybridged(area_target,g0)
%
% Make gradient waveform with desired area, constrained by specified edge values (symmetric)
%
% INPUTS:
%   area  - G/cm*sec
%   g0    - edge value (G/cm)
%   mxg   - G/cm
%   mxs   - G/cm/sec

dt  = 4e-6;       % sample duration (sec)
s = mxs * dt;     % max change in g per sample (G/cm)
if area_target < 0
	s = -s
end

% approximate waveform
ramp = g0;
g = [ramp ramp];
area = sum(g)*dt;
while area < area_target 
	ramp = [ramp (ramp(end)+s)];
	ramp(ramp>mxg) = mxg;
	ramp(ramp<-mxg) = -mxg;
	g = [ramp fliplr(ramp)];
	area = sum(g)*dt;
end

% scale to achieve exact area
gtmp = g-g(1)*ones(size(g));
areatmp = sum(gtmp)*dt;
area_targettmp = area_target - g(1)*length(g)*dt;
gtmp = gtmp/areatmp*area_targettmp;
g = gtmp + g(1)*ones(size(g));

return;
