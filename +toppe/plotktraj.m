function plotktraj(fname, varargin)
% function plotktraj(fname, varargin)
%
% Plot TOPPE .mod file.
%
% function plotktraj
% % Plots k space trajectory of a readout module
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
% Inputs:
%  fname     Readout .mod file name
% Options:
%  rainbowPlot    true/false

% Will plot the trajectory in 'readout.mod' file if no filename is input.

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


%% Set defaults 

% set filename to 'readout.mod' if not specified
if nargin == 0, fname = 'readout.mod'; end

% default to single color plot
arg.rainbowPlot = false;

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% Create var for plotting
[~,gx,gy,~] = toppe.readmod(fname);

% Create colormap vector to show the order in which the k-space points are
% acquired, starting at red, going to dark blue (in rainbow order)
if arg.rainbowPlot == true
    jf = .1;  % "jet fraction" what fraction of the jet colormap to add and then cut 
    cmap = jet(size(gx,1)+round(jf*size(gx,1)));      % use jet since it looks like a rainbow
    cmap = cmap(1:size(gx,1),:);  % chop off dark red to visualize better
    cmap = cmap(end:-1:1,:);      % invert order to start at red
end

%% Plot k-space trajectory
%  (rainbow trajectory only works for single shot)
figure;
if size(gx,2) == 1 && arg.rainbowPlot == true
    waveform_str = 'Waveform #1';
    xtraj = cumtrapz(gx);
    ytraj = cumtrapz(gy);
    scatter(xtraj,ytraj,9,cmap,'filled','DisplayName',waveform_str)
else
    for ishot = 1:size(gx,2)
        waveform_str = sprintf('Waveform #%d',ishot);
        xtraj = cumtrapz(gx(:,ishot));
        ytraj = cumtrapz(gy(:,ishot));
        plot(xtraj,ytraj,'.','DisplayName',waveform_str); hold on;
    end
end
xlabel('kx'); ylabel('ky');
hold off; title('k-space trajectory');
legend('Location', 'Best')

return;
