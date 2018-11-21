function [im imsos] = ift3(D, varargin)
% Centered inverse 3DFT of a 3D/2D single- or multi-coil data matrix.
%
% function [im imsos] = ift3(D, varargin)
%
% Example usage:
%  >> size(D) = [64 64 20 16]; ims = ift3(D);                      % 16 receive coils; returns multi-coil 3D images
%  >> size(D) = [64 64 16];    ims = ift3(D, 'type', 2d);          % 16 receive coils; returns multi-coil 2D images
%  >> size(D) = [64 64 20];    ims = ift3(D, 'dozfft', false);    % does not do ift along kz; returns 3D image
%
% $Id: ift3.m,v 1.4 2018/10/26 21:46:41 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/ift3.m,v $

import toppe.*
import toppe.utils.*

nd = ndims(D);
if nd > 4
	error('ndims(D) > 4');
end

%% parse and check inputs
arg.type   = '3d';                     % default
arg.dozfft = true;                     % default
arg = vararg_pair(arg, varargin);      % fill in values from varargin

if strcmp(arg.type, '3d') & nd < 3
	error('ndims(D) must be >= 3 for 3d recon');
end

if strcmp(arg.type, '2d') & nd > 3
	error('ndims(D) for 2d input must be no larger than 3');
end

%% number of receive channels (coils)
if strcmp(arg.type, '2d')
	ncoils = size(D, 3);
else
	ncoils = size(D, 4);
end
		
%% do ift
for coil = 1:ncoils
	if strcmp(arg.type, '2d')
		im(:,:,coil) = sub_ift3(D(:,:,coil), false);
	else
		im(:,:,:,coil) = sub_ift3(D(:,:,:,coil), arg.dozfft);
	end
end

imsos = sqrt(sum(abs(im).^2, nd));

return;


function im = sub_ift3(D, do3dfft)

if do3dfft
	im = fftshift(ifftn(fftshift(D)));
else
	% don't do fft in 3rd dimension
	for k = 1:size(D,3)
		im(:,:,k) = fftshift(ifftn(fftshift(D(:,:,k))));
	end
end

return;

