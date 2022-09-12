function [mask, X, Y, Z] = roi2mask(roi, N, FOV)
% function [mask, X, Y, Z] = roi2mask(roi, N, FOV)
%
% Represent a rectangular ROI by a mask on grid defined
% by N and FOV.
%
% Inputs:
%   roi     struct containing rectangular ROI parameters. See toppe.getroi()
%   N       [1 3]   matrix size
%   FOV     [1 3]   field of view (cm)
%   shoROI  true/false   
%
% Output:
%   mask    [N(1) N(2) N(3)] logical mask
%   X/Y/Z   grid point locations (cm)

% internal calculations are in mm
FOV = FOV*10;   % mm

if nargin < 4
    showROI = true;
end

% Distance between center of edge points
% (max voxel distance from center of FOV is FOVc/2)
FOVc = FOV - FOV./N;

[X, Y, Z] = ndgrid(linspace(1,-1,N(1))*FOVc(1)/2, ...
                   linspace(1,-1,N(2))*FOVc(2)/2, ...
                   linspace(1,-1,N(3))*FOVc(3)/2);

% get voxel locations inside rectangle (before rotating and translating)
masktmp = ones(N);
masktmp(X > roi.w/2 | X < -roi.w/2) = 0;
masktmp(Y > roi.h/2 | Y < -roi.h/2) = 0;
masktmp(Z > roi.t/2 | Z < -roi.t/2) = 0;
X = X(logical(masktmp));  % column vector
Y = Y(logical(masktmp));
Z = Z(logical(masktmp));

% rotate
LOCS = roi.rotmat * [X Y Z]';  % roi.rotmat = 3x3 rotation matrix
X = LOCS(1,:);
Y = LOCS(2,:);
Z = LOCS(3,:);

% translate
X = X + roi.x;
Y = Y + roi.y;
Z = Z + roi.z;

% remove locations outside FOV
inds = X > FOVc(1)/2 | X < -FOVc(1)/2 ...
     | Y > FOVc(2)/2 | Y < -FOVc(2)/2 ...
     | Z > FOVc(3)/2 | Z < -FOVc(3)/2;
X(inds) = [];
Y(inds) = [];
Z(inds) = [];

% convert voxel locations (in mm) to matrix indeces
% Note negative sign due to universal coordinate system (UCS)
Xinds = N(1)/2 +  round(-X/FOV(1)*N(1));     
Yinds = N(2)/2 +  round(-Y/FOV(2)*N(2));    
Zinds = N(3)/2 +  round(-Z/FOV(3)*N(3));    

% create mask
mask = zeros(N);
for ii = 1:length(Xinds)
    mask(Xinds(ii), Yinds(ii), Zinds(ii)) = 1;
end

% fill in gaps due to rounding
se = strel('cube', 3);
mask = imdilate(mask, se);
mask = imerode(mask, se);

mask = logical(mask);

% convert to cm
X = X/10; 
Y = Y/10; 
Z = Z/10; 

return

