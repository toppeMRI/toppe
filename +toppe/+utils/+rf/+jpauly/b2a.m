% function [aca] = b2a(bc);
%
% This function takes a b polynomial, and returns the minimum phase, 
%   minimum power a polynomial
%
% Inputs:
%   bc - beta polynomial coefficients
%
% Outputs:
%   aca - minimum phase alpha polynomial
%
% Original Code from John Pauly's rf_tools package
% Modified (slightly) by Peder Larson, 12/13/2005
% (c) Board of Trustees, Leland Stanford Junior University

function [aca] = b2a(bc);

import toppe.utils.rf.jpauly.*

n = length(bc);
% calculate minimum phase alpha
bcp = bc;
bl = length(bc);
blp = bl*8;
bcp(blp) = 0;
bf = fft(bcp);
bfmax = max(abs(bf));
if bfmax>=1.0,                % PM can result in abs(beta)>1, not physical
  bf = bf/(1e-8 + bfmax);  %   we scale it so that abs(beta)<1 so that
end;                          %   alpha will be analytic
afa = mag2mp(sqrt(1-bf.*conj(bf)));
aca = fft(afa)/blp;
aca = aca(n:-1:1);
aca = reshape(aca, size(bc));

