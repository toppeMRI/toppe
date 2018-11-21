function rf = makebs(amp, varargin)
% Create Fermi pulse and write to bs.mod
%
% function rf = makebs(amp, varargin)
%
% Inputs:
%   amp           amplitude, in Gauss
% Options
%   system        (optional) struct specifying hardware system info, see systemspecs.m
%
% See also calckbs.m
%
% $Id: makebs.m,v 1.1 2018/11/02 14:24:52 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+rf/makebs.m,v $

import toppe.*
import toppe.utils.*
import toppe.utils.rf.*

% parse inputs
arg.system = systemspecs();         % default
arg.showSpectrum = false;           % default
arg = vararg_pair(arg, varargin);   % substitute varargin values as appropriate

% Fermi pulse 
wrf = 0;   % frequency offset (Hz)
[rf,T] = fermi(amp,wrf,4e-3,0.1e-3);
rf = [0; rf(:); 0];   % make sure RF waveform starts and ends at 0 to avoid issues on scanner

% write mod file
writemod('rf', rf, 'ofname', 'bs.mod', ...
         'desc', 'Fermi pulse for Bloch-Siegert b1+ mapping', ...
         'system', arg.system );

if arg.showSpectrum
	spectrum = abs(fftshift(ifft(fftshift(rf))));
	k = (1/(2*T(end)))*(-length(T)/2:(length(T)/2)-1);
	k = k/1000;
	ii = find(abs(k)<max(2*abs(wrf),5));
	plot(k(ii),spectrum(ii));
	xlabel('kHz');
	title('Fermi pulse spectrum');
end

return;

