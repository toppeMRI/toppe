function plotktraj(fname)
% Plot TOPPE .mod file.
%
% function plotktraj
% % Plots k space trajectory of a readout module
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

[~,gx,gy,gz] = toppe.readmod(fname);

% Plot k-space trajectory
figure; 
for ishot = 1:size(gx,2)
    waveform_str = sprintf('Waveform #%d',ishot);
    xtraj = cumtrapz(gx(:,ishot));
    ytraj = cumtrapz(gy(:,ishot));
    plot(xtraj,ytraj,'.','DisplayName',waveform_str); hold on;
end
xlabel('kx'); ylabel('ky');
hold off; title('k-space trajectory');
legend('Location', 'Best')

return;