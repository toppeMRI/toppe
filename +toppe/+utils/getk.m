function [kxo, kxe] = getk(sysGE, readoutFile, nfid)
% Get kspace sample locations (assumes dwell = 4us)
% 
% Inputs:
%  sysGE
%  readoutFile    .mod file name containing ADC window
%  nfid           number of acquired samples in ADC window
%
% Outputs:
%  kxo            [1 nfid], odd echo k-space, cycles/cm
%  kxe            [1 nfid], even echo k-space, cycles/cm

[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readoutFile);
kx = sysGE.raster*1e-6*sysGE.gamma*1e-4*cumsum(gx);  % cycles/cm
kx = kx - kx(end)/2;
kx = kx((hdr.npre+1):(hdr.npre+nfid));
del = -0.0;
kxo = interp1(1:nfid, kx, (1:nfid) - 0.5 - del, 'linear', 'extrap');
kxe = interp1(1:nfid, kx, (1:nfid) + 0.5 + del, 'linear', 'extrap');
kxe = fliplr(kxe);
