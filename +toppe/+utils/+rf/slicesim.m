  function [m] = slicesim(m0, rf, gz, dt, Z, T1, T2, doDisplay)
% function [m] = slicesim(m0, rf, gz, dt, Z, T1, T2, [doDisplay])
%
% Simulate slice-selective rf pulse and display slice profile.
%
% Inputs:
%	 m0   [1 3]        initial magnetization (e.g., [0 0 1]);
%   rf   [nstep 1]    complex rf waveform (Gauss)
%   gz   [nstep 1]    gradient waveform (Gauss/cm)
%   dt   [1 1]        raster (sample) time (ms)
%   Z    [nz 1]       Spatial positions (cm) at which to place spins.
%   T1   [1 1]        ms
%   T2   [1 1]        ms
%   doDisplay  bool   plot slice profile? Default: true
%
% Output:
%   m   [length(Z) 1]   transverse magnetization, normalized to m0 (complex)
%
% Example:
%  >> m0 = [0 0 1];   % initial magnetization along mz
%  >> [rf,~,~,gz] = toppe.readmod('tipdown.mod');
%  >> Z = linspace(-5,5,100);       % cm
%  >> dt = 4e-3;                    % ms
%  >> toppe.utils.rf.slicesim(m0,rf,gz,dt,Z,T1,T2);

if strcmp(m0, 'test')
	sub_test();
	return;
end


if ~exist('doDisplay', 'var')
	doDisplay = true;
end

% work with column vectors
rf = rf(:);
gz = gz(:);
Z = Z(:);

dz = Z(2)-Z(1);   % distance between grid points

% effective field (z component)
Bz = gz*Z'*1e-4;     % [nstep nz], Tesla

nstep = length(rf);

for ii = 1:length(Z)
	Beff = [real(rf)*1e-4 imag(rf)*1e-4 Bz(:,ii)];                % [nstep 3], Tesla
	mtmp = toppe.utils.rf.blochsim(m0, Beff, T1, T2, dt, nstep);
	m(ii) = mtmp(end,1) + 1i*mtmp(end,2);
end

if doDisplay
T = dt*(1:nstep);
subplot(131); hold off; plot(T,abs(rf),'b'); hold on; plot(T,gz/20,'g'); legend('rf', 'gz/20'); xlabel('time (ms)');
subplot(132); plot(Z,abs(m)); xlabel('distance (cm)'); ylabel('abs(m)');
subplot(133); plot(Z,angle(m)); xlabel('distance (cm)'); ylabel('angle(m)');
end

return


function sub_test()

% design slice-selective RF pulse
flip = 90;        % degrees
slthick = 1.0;    % cm
tbw = 8; 
dur = 3;          % ms
ncycles = 4;
[rf, gz] = toppe.utils.rf.makeslr(flip, slthick, tbw, dur, ncycles, ...
	'type', 'ex', ...  % 90 degree excitation pulse (default is 'st' -- small tip)
	'ofname', 'test.mod');

% simulate slice profile
m0 = [0 0 1]; % initial magnetization vector
dt = 4e-3;    % ms
Z = linspace(-1*slthick, 1*slthick, 100);  % spatial positions at which to evaluate magnetization
T1 = 1000; T2 = 100;  % msec
toppe.utils.rf.slicesim(m0, rf, gz, dt, Z, T1, T2);

return;
