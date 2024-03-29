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

import toppe.*
import toppe.utils.*
import toppe.utils.rf.*

% parse inputs
arg.system = systemspecs();         % default
arg.showSpectrum = false;           % default
arg.ofname = 'bs.mod';
arg = vararg_pair(arg, varargin);   % substitute varargin values as appropriate

% Fermi pulse 
wrf = 0;   % frequency offset (Hz)
[rf,T] = fermi(amp,wrf,4e-3,0.1e-3);
rf = [0; rf(:); 0];   % make sure RF waveform starts and ends at 0 to avoid issues on scanner

rf = makeGElength(rf);

% write mod file
writemod(arg.system, ...
        'rf', rf, 'ofname', arg.ofname, ...
        'desc', 'Fermi pulse for Bloch-Siegert b1+ mapping');

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

