function g = makeGElength(g)
% If length not divisible by 2, pad with zeroes at end
%
% function g = makeGElength(g)
%
% Input:
%   g    waveform array [nt ...]
%

[res, npulses] = size(g, 1);
g = [g; zeros(mod(res, 2), npulses)];

