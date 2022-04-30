function d = loadepi(pfile, echo, readoutFile)
% function d = loadbssfp4ge(pfile, echo, readoutFile)
%
% Load raw data Pfile acquired with scan created by makeepi.m
%
% Inputs:
%   pfile         string
%   echo          int
%   readoutFile   string
%
% Output:
%   d    [nx*decimation ny nz nCoils]

% load data
din = toppe.utils.loadpfile(pfile, echo);  % [ndat nCoils nslices 1 nviews]
din = flipdim(din, 1); % !! Data seems to always be flipped along readout (do check though)
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

