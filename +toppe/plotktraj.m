function plotktraj(fname, varargin)
% Plot TOPPE .mod file.
%
% function plotktraj
% % Plots k space trajectory of a readout module
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.

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

%%% convert to physical units. Assumes gradients are in Gauss/cm, gamma is
%%% gamma/2pi in Hz/Gauss, and raters is in s. If gradients are in other
%%% units would need to change these values.
sys = toppe.systemspecs();
gamma = sys.gamma;   % gamma/2pi in Hz/Gauss
raster = 1e-6*sys.raster; % sampling time in s

% calculate and plot ktraj
figure;
if size(gx,2) == 1 && arg.rainbowPlot == true
    kxtraj = (gamma * cumtrapz(gx) * raster) * 0.1; % kxtraj in 1/mm
    kytraj = (gamma * cumtrapz(gy) * raster) * 0.1; % kytraj in 1/mm
    scatter(kxtraj,kytraj,9,cmap,'filled')
    axis([1.1*min(kxtraj), 1.1*max(kxtraj), 1.1*min(kytraj), 1.1*max(kytraj)])
    axis equal
else
    if arg.rainbowPlot == true
        for ishot = 1:size(gx,2)
            kxtraj = (gamma * cumtrapz(gx(:,ishot)) * raster) * 0.1; % kxtraj in 1/mm
            kytraj = (gamma * cumtrapz(gy(:,ishot)) * raster) * 0.1; % kytraj in 1/mm
            scatter(kxtraj,kytraj,9,cmap,'filled')
            hold on;
        end
    else
        for ishot = 1:size(gx,2)
            waveform_str = sprintf('Waveform #%d',ishot);
            kxtraj = (gamma * cumtrapz(gx(:,ishot)) * raster) * 0.1; % kxtraj in 1/mm
            kytraj = (gamma * cumtrapz(gy(:,ishot)) * raster) * 0.1; % kytraj in 1/mm
            plot(kxtraj,kytraj,'.','DisplayName',waveform_str); hold on;
            hold on;
        end
        legend('Location', 'Best')
    end
    kxtraj_all = (gamma * cumtrapz(gx) * raster) * 0.1;
    kytraj_all = (gamma * cumtrapz(gy) * raster) * 0.1;
    axis([1.2*min(kxtraj_all(:)), 1.2*max(kxtraj_all(:)), 1.2*min(kytraj_all(:)), 1.2*max(kytraj_all(:))])
    axis equal
    
end
xlabel('kx (mm^{-1})'); ylabel('ky (mm^{-1})');
hold off; title('k-space trajectory');


return;
