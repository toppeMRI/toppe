function [kspace, acq_indices, grad_amps] = calckspace(fname_readmod, varargin)

% [kspace, acq_indices, grad_amps] = calckspace(fname, varargin)

% Optional inputs:
%   fname      filename of the readout module

%
% Outputs:
%   kspace       [Nk n_unique_shots 3] kspace samples in units 1/mm
%   acq_indices  [total_shots 4]       indices of slice, echo, view & rot at
%                                      each shot
%   grad_amps    [total_shots 3]       gx, gy, gz amplitude at each shot

%%% NOTE! the echo index in toppe starts at 0, here it starts at 1 to make
%%% it MATLAB indexing compatible

% Melissa Haskell, mhask@umich.edu, June 2020

%% parse inputs and set defaults
arg.loopFile        = 'scanloop.txt';
arg.system          = toppe.systemspecs();  % Accept default timing (includes EPIC-related time gaps)
arg = toppe.utils.vararg_pair(arg, varargin);

if nargin < 1, fname_readmod = 'readout.mod'; end

%% Set constants
max_pg_iamp = 2^15-2;

%% Convert to physical units. Assumes:
%   -gradients are in Gauss/cm, gamma is
%   -gamma/2pi in Hz/Gauss
%   -raster is in s. 
% If gradients are in other units would need to change these values.

sys = arg.system;
if ~strcmp(sys.gradUnit,'Gauss/cm')
    disp(['kspace calculated assuming gradients in units of Gauss/cm,',...
        ' but system units are different. k-space values may be inaccurate.'])
end
gamma = sys.gamma;   % gamma/2pi in Hz/Gauss
raster = sys.raster; % sampling time in s

%% Calculate and return kspace for readout module

% gradients in Gauss/cm
[~,gx,gy,gz] = toppe.readmod(fname_readmod);
% n_unique_shots = size(gx,2);
% if size(gx,2) ~= size(gy,2), warn('check readout.mod sizes'); end

% Based on equation (5.36) in Nishimura Principles of MRI book
kxtraj = (gamma * cumtrapz(gx) * raster) * 0.1; % kxtraj in 1/mm
kytraj = (gamma * cumtrapz(gy) * raster) * 0.1; % kytraj in 1/mm
kztraj = (gamma * cumtrapz(gz) * raster) * 0.1; % kztraj in 1/mm

kspace = cat(3, kxtraj,kytraj, kztraj);

%% read scanloop and get info needed for k-space recon
loopArr = toppe.tryread(@toppe.readloop, arg.loopFile);

% extract loop rows where data is saved to P-file
loopArr_dabslice = loopArr(:,7);
loopArr_dabmode = loopArr(:,10);
loopArr_daq = loopArr(loopArr_dabslice & loopArr_dabmode,:);

% get gradient amplitude scales
grad_amps = loopArr_daq(:,4:6) / max_pg_iamp;

% get slice, echo, shot, and rotation ordering
acq_indices = loopArr_daq(:,[7:9,11]);

% start view indexing at 1
acq_indices(:,2) = acq_indices(:,2) + 1; 

% ntotal_shots = size(acq_indices,1);
% 
% kspace = zeros(size(kxtraj,1), 3, ntotal_shots);
% for ii = 1:ntotal_shots
% 
%     shot_rot = loopArr_daq(ii,11);
%     shot_num = loopArr_daq(ii,9);
%     kxy_tmp = kxtraj(:,shot_num) + 1i*kytraj(:,shot_num);
%     kxy_shot = kxy_tmp * exp(-1i*shot_rot);
%     kspace(:,:,ii) = [real(kxy_shot), imag(kxy_shot), kztraj(:,shot_num)];
% end


end


