function [gxsame, gysame, gzsame] = makesamelength(gx,gy,gz)
% Makes waveforms the same length by padding zeros to shorter ones
% function g = makesamelength(gx,gy,gz) (gz optional)
% In: gx,gy,gz [nsamp nwaveforms] 

if nargin < 2
    error('Must specify at least two waveforms');
elseif nargin == 2
    gz = zeros(size(gx));
end

g_length(1) = size(gx,1);
g_length(2) = size(gy,1);
g_length(3) = size(gz,1);
g_maxlength = max(g_length);

nwaveforms(1) = size(gx,2);
nwaveforms(2) = size(gy,2);
nwaveforms(3) = size(gz,2);

if ~all(nwaveforms == nwaveforms(1))
    error('Inputs must have the same number of waveforms (dim 2)');
else
    nwaveforms = nwaveforms(1);
end

gxsame = [gx; zeros(g_maxlength-g_length(1),nwaveforms)];
gysame = [gy; zeros(g_maxlength-g_length(2),nwaveforms)];
gzsame = [gz; zeros(g_maxlength-g_length(3),nwaveforms)];
return
