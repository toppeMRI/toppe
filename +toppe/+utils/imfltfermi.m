function Dout = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
% Apply low-pass 2D Fermi filter (circular or square) to each slice in a 3D image volume
%
% function Dout = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
%
% D                   2D/3D/nD image volume. Can also be kspace, but in that case set Dtype = 'kspace'.
%                     If dimension is > 3, it is assumed that D contains multiple 3D volumes.
% fltradius           radius of filter, in pixels
% transitionwidth     Fermi filter roll-off (75%-25% width)
% fltgeom             'circ' (default) or 'rect'
% Dtype               'image' (default) or 'kspace'

import toppe.utils.*

nims = size(D,4);
sz = size(D);

D = reshape(D, size(D,1), size(D,2), size(D,3), []);

if ~exist('fltgeom', 'var')
	fltgeom = 'circ';
end
if ~exist('Dtype', 'var')
	Dtype = 'image';
end

flt = fermi2d(size(D,1), fltradius, transitionwidth, fltgeom);

for ii = 1:nims
	for iz = 1:size(D,3)
		if strcmp(Dtype, 'image')
			d = cfftn(D(:,:,iz,ii), 'forward');
		else
			d = D(:,:,iz,ii);
		end
		Dout(:,:,iz,ii) = bsxfun(@times, d, flt);
	end
end

if strcmp(Dtype, 'image')
	for ii = 1:nims
		Dout(:,:,:,ii) = cfftn(Dout(:,:,:,ii), 'inverse');
	end
end

Dout = reshape(Dout, sz);

return;
