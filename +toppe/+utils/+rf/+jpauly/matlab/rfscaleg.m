function rfs = rfscaleg(rf, pulseduration)
% RFSCALEG Converts a normalized RF waveform to Gauss
%   RF_G = RFSCALEG(RF, PW)
% Inputs:
%   RF - scaled to sum to the desired flip angle in radians
%   PW - Pulse duration in ms
% Outputs:
%   RF_G - rf pulse scaled to Gauss
%
% Peder Larson 10/6/04

GAMMA = 4257; % Hz/G

rfs = rf*length(rf)/(2*pi*GAMMA*pulseduration*1e-3);