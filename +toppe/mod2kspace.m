function [kspace] = mod2kspace(fname, varargin)
% load a TOPPE .mod file and create kspace based on system specs.
%
% function mod2kspace
%
% Inputs:
%   fname - filename of readout module. Default is 'readout.mod'
%   varargin - variable arguments
%
% Outputs:
%   kspace - 2d or 3d k-space trajectory in units 1/mm
%
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.

% Will calculate the kspace in 'readout.mod' file if no filename is input.

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


%% Set defaults and parse inputs

% set filename to 'readout.mod' if not specified
if nargin == 0, fname = 'readout.mod'; end

% Substitute varargin values as appropriate (unused as of June 12, 2020)
if nargin > 1
    arg = toppe.utils.vararg_pair(arg, varargin);
end

%% Load gradients from mod file.
[~,gx,gy,gz] = toppe.readmod(fname);


%% Convert to physical units. Assumes:
%   -gradients are in Gauss/cm, gamma is
%   -gamma/2pi in Hz/Gauss
%   -raster is in s. 
% If gradients are in other units would need to change these values.

sys = toppe.systemspecs();
if ~strcmp(sys.gradUnit,'Gauss/cm')
    disp(['kspace calculated assuming gradients in units of Gauss/cm,',...
        ' but system units are different. k-space values may be inaccurate.'])
end
gamma = sys.gamma;   % gamma/2pi in Hz/Gauss
raster = sys.raster; % sampling time in s

%% Calculate and return kspace

% Based on equation (5.36) in Nishimura Principles of MRI book
kxtraj = (gamma * cumtrapz(gx) * raster) * 0.1; % kxtraj in 1/mm
kytraj = (gamma * cumtrapz(gy) * raster) * 0.1; % kytraj in 1/mm
kztraj = (gamma * cumtrapz(gz) * raster) * 0.1; % kztraj in 1/mm

% only return 2 columns if no gz gradients played
if nnz(kztraj) == 0
    kspace = [kxtraj(:), kytraj(:)];
else
    kspace = [kxtraj(:), kytraj(:), kztraj(:)];
end

return;