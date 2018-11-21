function sub_plotmod(fname)
% Plot TOPPE .mod file.
%
% function sub_plotmod(fname)
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.

% TOPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU Library General Public License as published by
% the Free Software Foundation version 2.0 of the License.
%
% TOPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public License
% along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
% 
% (c) 2016 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

[b1,gx,gy,gz] = toppe.readmod(fname);
rho = abs(b1);
theta = angle(b1);

nt = size(b1,1);
dt = 4e-3;  % ms
T = linspace(dt/2,nt*dt-dt/2,nt);

figure
% rf
subplot(337); sub_plot(T, rho);   ylabel('abs(rf) G');
xlabel('time (msec)');
subplot(338); sub_plot(T, theta); ylabel('angle(rf) rad');
xlabel('time (msec)');

% gradient waveform
subplot(331); sub_plot(T,gx);    ylabel('gx G/cm');
subplot(332); sub_plot(T,gy);    ylabel('gy G/cm');
subplot(333); sub_plot(T,gz);    ylabel('gz G/cm');

% gradient slew
subplot(334); sub_plot(T(2:end), diff(gx)/dt);  ylabel('gx slew    G/cm/ms');
subplot(335); sub_plot(T(2:end), diff(gy)/dt);  ylabel('gy slew    G/cm/ms');
subplot(336); sub_plot(T(2:end), diff(gz)/dt);  ylabel('gz slew    G/cm/ms');
xlabel('time (msec)');

return;

function sub_plot(T,wav,label)

h = line(T,wav);
set(h, 'LineWidth', 1.0);

return;

