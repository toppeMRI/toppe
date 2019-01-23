function rf = b2rf(bc);
% function rf = b2rf(bc);
%
% This function takes a beta polynomial and returns an RF pulse
% A minimum phase alpha polynomial is intermediately computed,
% followed by the SLR transform.
% 
% Works for complex RF pulses
% Requires b2a.m and ab2rf.m
%
% Inputs:
%   bc - beta polynomial coefficients
%
% Outputs:
%   rf - RF pulse
%
% Original Code from John Pauly's rf_tools package
% Modified (slightly) by Peder Larson, 12/13/2005
% (c) Board of Trustees, Leland Stanford Junior University

% Compute minimum phase alpha
ac = b2a(bc);
rf = ab2rf(ac, bc);