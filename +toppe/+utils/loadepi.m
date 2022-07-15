function [d, kx, ky, kz, ti] = loadepi(pfile, echo, readoutFile, varargin)
% function [d, kx, ky, kz, ti] = loadepi(pfile, echo, readoutFile, varargin)
%
% Load raw data Pfile acquired with scan created by makeepi.m
%
% Inputs:
%   pfile         string
%   echo          int
%   readoutFile   string
%
% Keyword-argument input options:
%   flipFirstDim  true/false   Flip data along the FID dimension? Default: true.
%   zfov          [1]  cm. Default: 24
%
% Outputs:
%   d    [nx*decimation ny nz nCoils]
%        In TOPPE, samples are acquired every 4us.
%        The 'effective'/design dwell time is decimation*4us, i.e.,
%        the image needs to be cropped (in x) if decimation > 1.
%   kx/ky   [nx*decimation ny nz]  cycles/cm
%   ti   [nx*decimation ny nz]  sample times (sec)

arg.flipFID = true;
arg.zfov = 24;   

arg = vararg_pair(arg, varargin);      % requires MIRT

% load data
din = toppe.utils.loadpfile(pfile, echo);  % [ndat nCoils nslices 1 nviews]
if arg.flipFID
    din = flipdim(din, 1); 
end
din = squeeze(din); % [nFID nCoils nz nShots]

[nFID nCoils nz nShots] = size(din);

% get readout parameters
% see makeepi.m
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutFile);
nx = paramsint16(1); 
ny = paramsint16(2); 
Ry = paramsint16(4);
nSampPerEcho = paramsint16(5) + paramsint16(6);
decimation = paramsint16(7);
nRamp = paramsint16(8);
etl = ny/nShots/Ry;  % echo train length

% kspace
% apply 1/2 sample shift which seems about right for EPI
raster = 4e-6;      % sec
gamma = 4.2576e3;   % Hz/T
kx = raster*gamma*cumsum(gx,1);  % cycles/cm
ky = raster*gamma*cumsum(gy,1);
nk = size(kx,1);
kx = interp1(1:nk, kx, (1:nk) - 0.5, 'linear', 'extrap');
ky = interp1(1:nk, ky, (1:nk) - 0.5, 'linear', 'extrap');
kx = kx - max(kx)/2; % center at k=0 (assume the prephaser does this)
ky = ky - max(ky)/2;
kz = zeros([nx*decimation ny nz]);
dkz = 1/zfov;  % cycles/cm  (kspace spacing)
for iz = 1:nz
    kz(:,:,iz) = (iz-nz/2-0.5)*nz*dkz;
end

% sample times
ti = raster*[[1:nk] - 0.5];  % sec

% sort data into [nx*decimation ny nz nCoils] matrix
fprintf('Reshaping data matrix... ');
d = zeros(nx*decimation, nCoils, nz, ny);
for ish = 1:nShots
    for e = 1:etl  % loop over echoes
        iy = ish + (e-1)*nShots;
        iStart = (e-1)*nSampPerEcho + nRamp + 1;
        iStop = iStart + nx*decimation - 1;
        d(:, :, :, iy) = din(iStart:iStop, :, :, ish);   
    end
end
d = permute(d, [1 4 3 2]);
fprintf('done\n');

return

