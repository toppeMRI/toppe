function rf = makehardpulse(flip, dur, raster)
% function rf = makehardpulse(flip, dur, raster = 4e-6)
%
% Output waveform is padded with zero at beginning and end.
%
% flip      degrees
% dur       sec
% raster    sec (default: 4e-6)
%
 
if nargin < 3
	raster = 4e-6;
end

gamma = 4.2576e+03;   % Hz/Gauss
n = 2*round(dur/raster/2);
dur = n*raster;
flip = flip/180*pi;    % radians
rfamp = flip / (2*pi*gamma*dur);    % Gauss
rf = [zeros(1,1); ones(n,1)*rfamp; zeros(1,1)];

return
