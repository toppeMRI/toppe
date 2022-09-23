function [X, Y, Z] = grid2xyz(N, FOV)
% function [X, Y, Z] = grid2xyz(N, FOV)
%
% Return spatial locations consistent with the 
% 'universal coordinate system' (UCS) that is used by various
% utilities such as the Java SlicePlanner and the 
% B0 shimming tool.
% This is a right-handed system with cm units:
% Patient R = +x; A = +y; S = +z; with (0,0,0) at scanner iso-center

FOVc = FOV - FOV./N;   % Distance between edge grid points

[X, Y, Z] = ndgrid(linspace(1,-1,N(1))*FOVc(1)/2, ...
                   linspace(1,-1,N(2))*FOVc(2)/2, ...
                   linspace(1,-1,N(3))*FOVc(3)/2);

