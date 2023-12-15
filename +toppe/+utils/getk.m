function [kxo, kxe] = getk(sysGE, readoutFile, nfid, del)
% Get kspace sample locations from .mod file
% 
% Inputs:
%  sysGE          see toppe.systemspecs()
%  readoutFile    .mod file name containing ADC window
%  nfid           number of acquired samples in ADC window
%
% Option:
%  del            offset by this many samples (default: 0)
%
% Outputs:
%  kxo            [1 nfid], odd echo k-space, cycles/cm
%  kxe            [1 nfid], even echo k-space, cycles/cm

if nargin < 4
    del = 0.0;
end

[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutFile);
kx = sysGE.raster*1e-6*sysGE.gamma*1e-4*cumsum(gx);  % cycles/cm
kx = kx - kx(end)/2;
kx = kx((hdr.npre+1):(hdr.npre+nfid));
kxo = interp1(1:nfid, kx, (1:nfid) - 0.5 - del, 'linear', 'extrap');
kxe = interp1(1:nfid, kx, (1:nfid) + 0.5 + del, 'linear', 'extrap');
kxe = fliplr(kxe);
