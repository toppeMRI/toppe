function [ims rss] = ift3(D, varargin)
% Centered inverse 3DFT of a 3D/2D single- or multi-coil data matrix.
%
% function [ims rss] = ift3(D, varargin)
%
% Input
%   D      [nx ny (nz) nCoils]
%
% Optional inputs
%   decimation   int    effective dwell time is decimation*4us. Default: 1
%   dozfft       true/false   Default: true
%   type         '2d'/'3d'  Specifies whether a 3D input matrix is treated
%                as 3D image or 2D multicoil images. Default: '3d'
% Example usage:
%  >> size(D) = [64 64 20 16]; ims = ift3(D);                      % 16 receive coils; returns multi-coil 3D images
%  >> size(D) = [64 64 20];    ims = ift3(D, 'type', 2d);          % 20 receive coils; returns multi-coil 2D images
%  >> size(D) = [64 64 20];    ims = ift3(D, 'dozfft', false);     % does not do ift along kz; returns 3D image
%

import toppe.*
import toppe.utils.*

% parse inputs and dimensions
arg.decimation = 1; 
arg.dozfft = true; 
arg.type   = '3d';
arg = vararg_pair(arg, varargin);      % fill in values from varargin

nd = ndims(D);

if strcmp(arg.type, '3d') & nd < 3
	error('ndims(D) must be >= 3 for 3d recon');
end

if strcmp(arg.type, '2d') & nd > 3
	error('ndims(D) for 2d input must be no larger than 3');
end

if nd > 4 | nd < 2
	error('D must be 2D or 3D single/multi-coil matrix');
end

nx = size(D,1);
ny = size(D,2);

if strcmp(arg.type, '2d')
    nCoils = size(D, 3);
    nz = [];
else
    nCoils = size(D, 4);
    nz = size(D, 3);
end

% do ift
for coil = 1:nCoils
	if strcmp(arg.type, '2d')
		ims(:,:,coil) = sub_ift3(D(:,:,coil), false);
	else
		ims(:,:,:,coil) = sub_ift3(D(:,:,:,coil), arg.dozfft);
	end
end

% Trim according to oversampling in x
nx = size(ims,1);
ims = ims(end/2+((-nx/arg.decimation/2):(nx/arg.decimation/2-1))+1,:); 
ims = reshape(ims, [nx/arg.decimation ny nz nCoils]);

% root sum of squares coil combination
if nCoils > 1
    rss = sqrt(sum(abs(ims).^2, nd));
else
    rss = abs(ims);
end

return;


function ims = sub_ift3(D, do3dfft)

if do3dfft
	ims = fftshift(ifftn(fftshift(D)));
else
	% don't do fft in 3rd dimension
	for k = 1:size(D,3)
		ims(:,:,k) = fftshift(ifftn(fftshift(D(:,:,k))));
	end
end

return;

