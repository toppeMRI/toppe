function [im imsos] = ft3(D, varargin)
% Centered forward 3DFT of a 3D/2D single- or multi-coil data matrix.
%
% function [im imsos] = ft3(D)

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
		
%% do ft
for coil = 1:ncoils
	if strcmp(arg.type, '2d')
		im(:,:,coil) = sub_ft3(D(:,:,coil), false);
	else
		im(:,:,:,coil) = sub_ft3(D(:,:,:,coil), arg.dozfft);
	end
end

imsos = sqrt(sum(abs(im).^2, nd));

return;


function im = sub_ft3(D, do3dfft)

if do3dfft
	im = fftshift(fftn(fftshift(D)));
else
	% don't do fft in 3rd dimension
	for k = 1:size(D,3)
		im(:,:,k) = fftshift(fftn(fftshift(D(:,:,k))));
	end
end

return;

