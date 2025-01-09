function R = angleaxis2rotmat(alpha, u)
% function R = angleaxis2rotmat(alpha, u)
%
% Inputs
%  alpha   [1]              Rotation angle (radians)
%  u       [3 1] or [1 3]   Vector lying along the axis of rotation
%
% Output
%  R       [3 3] rotation matrix
%
% Useful for, e.g., applying "in-plane" rotations for multi-shot spiral imaging.
%
% From https://www.mathworks.com/matlabcentral/fileexchange/66446-rotation-matrix?s_tid=mwa_osa_a
% Accessed 22-Sep-2019

% RotMatrix - N-dimensional Rotation matrix
% R = RotMatrix(alpha, u, v)
% INPUT:
%   alpha: Angle of rotation in radians, counter-clockwise direction.
%   u, v:  Ignored for the 2D case.
%          For the 3D case, u is the vector to rotate around.
%          For the N-D case, there is no unique axis of rotation anymore, so 2
%          orthonormal vectors u and v are used to define the (N-1) dimensional
%          hyperplane to rotate in.
%          u and v are normalized automatically and in the N-D case it is cared
%          for u and v being orthogonal.
% OUTPUT:
%   R:     Rotation matrix.
%          If the u (and/or v) is zero, or u and v are collinear, The rotation
%          matrix contains NaNs.
%
% REFERENCES:
% analyticphysics.com/Higher%20Dimensions/Rotations%20in%20Higher%20Dimensions.htm
% en.wikipedia.org/wiki/Rotation_matrix
% application.wiley-vch.de/books/sample/3527406204_c01.pdf
%
% Tested: Matlab 7.7, 7.8, 7.13, 9.1, WinXP/32, Win7/64
% Author: Jan Simon, Heidelberg, (C) 2018 matlab.2010(a)n(MINUS)simon.de

% $JRev: R-b V:001 Sum:GwB8LUFcZ+7i Date:10-Mar-2018 19:26:01 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GL3D\RotMatrix.m $
% History:
% 001: 10-Mar-2018 14:31, First version.

% Initialize: ==================================================================
% Global Interface: ------------------------------------------------------------
% Initial values: --------------------------------------------------------------
% Program Interface: -----------------------------------------------------------
if numel(alpha) ~= 1
   error('JSimon:RotMatrrix:BadInput1', ...
      'Angle of rotation must be a scalar.');
end

% User Interface: --------------------------------------------------------------
% Do the work: =================================================================
s = sin(alpha);
c = cos(alpha);

% Different algorithms for 2, 3 and N dimensions:
switch nargin
   case 1
      % 2D rotation matrix:
      R = [c, -s;  s, c];
      
   case 2
      if numel(u) ~= 3
         error('JSimon:RotMatrrix:BadAxis2D', ...
            '3D: Rotation axis must have 3 elements.');
      end
      
      % Normalized vector:
      u = u(:);
      u = u ./ sqrt(u.' * u);
            
      % 3D rotation matrix:
      x  = u(1);
      y  = u(2);
      z  = u(3);
      mc = 1 - c;
      R  = [c + x * x * mc,      x * y * mc - z * s,   x * z * mc + y * s; ...
            x * y * mc + z * s,  c + y * y * mc,       y * z * mc - x * s; ...
            x * z * mc - y * s,  y * z * mc + x * s,   c + z * z .* mc];
         
      % Alternative 1 (about 60 times slower):
      % R = expm([0, -z,  y; ...
      %           z,  0, -x; ...
      %          -y,  x,  0]  * alpha);
      
      % Alternative 2:
      % R = [ 0, -z, y; ...
      %       z, 0, -x; ...
      %      -y, x,  0] * s + (eye(3) - u * u.') * c + u * u.';
              
   case 3
      n = numel(u);
      if n ~= numel(v)
         error('JSimon:RotMatrrix:BadAxes3D', ...
            'ND: Axes to define plane of rotation must have the same size.');
      end
      
      % Normalized vectors:
      u = u(:);
      u = u ./ sqrt(u.' * u);
      
      % Care for v being orthogonal to u:
      v = v(:);
      v = v - (u.' * v) * u;
      v = v ./ sqrt(v.' * v);
      
      % Rodrigues' rotation formula:
      R = eye(n) + ...
         (v * u.' - u * v.') * s + ...
         (u * u.' + v * v.') * (c - 1);
      
   otherwise
      error('JSimon:RotMatrrix:BadNInput', ...
            '1 to 3 inputs required.');
end

end
