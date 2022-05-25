function [kbs th] = calckbs(rf, freq, dt)
% function [kbs th] = calckbs(rf, freq, dt)
%
% Calculate 'K_BS' for Bloch-Siegert mapping. See Sacolick et al MRM 2010 and example below.
%
% Inputs:
%   rf            rf waveform, Gauss
%   freq          off-resonance transmit frequency (Hz)
%   dt            rf sample/raster time (sec)
%
% Example:
%  rf = makebs(0.05);
%  dt = 4e-6;  % sec
%  kbs = calckbs(rf, 4000, dt);
%  % acquire data at off-resonance frequencies +4000 Hz and -4000 Hz
%  pc = phasecontrast3d(imPlus4000, imMinus4000);
%  b1 = sqrt(pc/2/kbs);   % Measured b1 in Gauss. Note factor of 2!
%

% expected Bloch-Siegert shift from this pulse, for a given freq
gamma = 4.2576e3;         % Hz/Gauss
th = 2*pi*sum( (gamma*abs(rf)).^2/(2*freq) ) * dt; % / pi*180

% K_BS
kbs = th/(max(abs(rf(:)))^2);

return;

