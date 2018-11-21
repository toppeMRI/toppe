function [rf,T] = fermi(A,wrf,t0,a)
% See MRI handbook, Eq. (4.14)
%
% $Id: fermi.m,v 1.1 2018/11/02 14:24:52 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+rf/fermi.m,v $

dt = 4e-6;   % sec
T = [-(t0+5*a):dt:(t0+5*a)];
rf = A*exp(i*2*pi*wrf*T) ./ ( 1 + exp((abs(T)-t0)/a) );

return;
