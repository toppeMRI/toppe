% ab2rf - take two polynomials for alpha and beta, and return an
%   rf waveform that would generate them.  
%
% COMPLEX RF VERSION
%
% rf = ab2rf(a, b)
%   a, b - polynomials for alpha and beta
%   rf - rf waveform that produces alpha and beta under the hard pulse
%     approximation.
%
% Courtesy John Pauly, 2004
% (c) Board of Trustees, Leland Stanford Junior University

function [rf] = ab2rf(ac,bc)

n = length(ac);
j = sqrt(-1);
for i=n:-1:1,
  c(i) = sqrt(1/(1+abs(bc(i)/ac(i))^2));
  s(i) = conj(c(i)*bc(i)/ac(i));
  theta = atan2(abs(s(i)),c(i));
  psi = angle(s(i));
  rf(i) = 2*(theta*cos(psi)+j*theta*sin(psi));
  acn = c(i)*ac + s(i)*bc;
  bcn = -conj(s(i))*ac + c(i)*bc;
  ac = acn(2:i);
  bc = bcn(1:i-1);
end;

