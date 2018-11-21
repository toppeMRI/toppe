function [kbs th] = calckbs(rf, freq)
% Calculate 'K_BS' for Bloch-Siegert mapping. See Sacolick et al MRM 2010 and example below.
%
% Inputs:
%   rf            rf waveform, Gauss
%   freq          off-resonance transmit frequency (Hz)
%
% Example:
%  rf = makebs(0.05);
%  kbs = calckbs(rf, 4000);
%  % acquire data at off-resonance frequencies +4000 Hz and -4000 Hz
%  pc = phasecontrast3d(imPlus4000, imMinus4000);
%  b1 = sqrt(pc/2/kbs);   % Measured b1 in Gauss. Note factor of 2!
%
% $Id: calckbs.m,v 1.2 2018/11/02 14:30:29 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+rf/calckbs.m,v $

import toppe.*

sys = systemspecs();

% expected Bloch-Siegert shift from this pulse, for a given freq
freq = 4000;   % Hz
dt = sys.raster;            % (sec) duration of each rf waveform sample
gamma = sys.gamma;       % Hz/Gauss
th = 2*pi*sum( (gamma*abs(rf)).^2/(2*freq) ) * dt; % / pi*180

% K_BS
kbs = th/(max(abs(rf(:)))^2);

return;

