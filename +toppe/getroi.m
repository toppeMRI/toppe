function roi = getroi(fname, roiId)
% function roi = getroi(fname, roiId)
%
% Load slice prescription (ROI) from HDF5 file. See ROI.java.
%
% Inputs
%   fname    [string]         HDF5 file containing one (or more) ROIs. Created with the GUI.
%   roiId    [1 1] integer    ROI index, starting at 1. Default: 1
%
% Output
%   roi      Struct containing 3D ROI dimensions, center placement, and rotation. See also ROI.java.

if ~exist('roiId', 'var')
	roiId = 1;
end

% ROI dimensions
roi.w = h5read(fname, sprintf('/ROI%d/dimensions/width',     roiId));
roi.h = h5read(fname, sprintf('/ROI%d/dimensions/height',    roiId));
roi.t = h5read(fname, sprintf('/ROI%d/dimensions/thickness', roiId));

% ROI size limits
roi.wmin = h5read(fname, sprintf('/ROI%d/dimensions/minWidth',     roiId));
roi.hmin = h5read(fname, sprintf('/ROI%d/dimensions/minHeight',    roiId));
roi.tmin = h5read(fname, sprintf('/ROI%d/dimensions/minThickness', roiId));
roi.wmax = h5read(fname, sprintf('/ROI%d/dimensions/maxWidth',     roiId));
roi.hmax = h5read(fname, sprintf('/ROI%d/dimensions/maxHeight',    roiId));
roi.tmax = h5read(fname, sprintf('/ROI%d/dimensions/maxThickness', roiId));

% center location
roi.x = h5read(fname, sprintf('/ROI%d/center/x', roiId));
roi.y = h5read(fname, sprintf('/ROI%d/center/y', roiId));
roi.z = h5read(fname, sprintf('/ROI%d/center/z', roiId));

% rotation matrix
rotv = h5read(fname, sprintf('/ROI%d/rotmat', roiId));   % vectorized, row-major order
for ii = 1:3
	for jj = 1:3
		roi.rotmat(ii,jj) = rotv((ii-1)*3+jj);
	end
end

roi.scanPlaneToIsocenterDistance = h5read(fname, sprintf('/ROI%d/scanPlaneToIsocenterDistance', roiId));

%varargout{1} = 1;
%varargout{2} = 2;
%varargout{3} = 3;

