function Dflt = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
% Apply low-pass 2D Fermi filter (circular or square) to each slice in a 3D image volume
%
% function Dflt = imfltfermi(D, fltradius, transitionwidth, fltgeom, Dtype)
%
% D                   2D/3D/nD image volume. Can also be kspace, but in that case set Dtype = 'kspace'.
%                     If dimension is > 2, it is assumed that D contains one or more 3D volumes.
% fltradius           radius of filter, in pixels
% transitionwidth     Fermi filter roll-off (75%-25% width)
% fltgeom             'circ' (default) or 'rect'
% Dtype               'image' (default) or 'kspace'

if strcmp(D, 'test')
	Dflt = sub_test;
	return;
end

nims = size(D,4);
sz = size(D);

D = reshape(D, size(D,1), size(D,2), size(D,3), []);

if ~exist('fltgeom', 'var')
	fltgeom = 'circ';
end
if ~exist('Dtype', 'var')
	Dtype = 'image';
end

% k-space filter (2D)
flt = toppe.utils.fermi2d(size(D,1), fltradius, transitionwidth, fltgeom);

% apply filter
for ii = 1:nims
	for iz = 1:size(D,3)
		if strcmp(Dtype, 'image')
			d = toppe.utils.cfftn(D(:,:,iz,ii), 'forward');         % 2D kspace
			dflt = bsxfun(@times, d, flt);                          % filtered
			Dflt(:,:,iz,ii) = toppe.utils.cfftn(dflt, 'inverse');   % filtered 2D image
		else
			d = D(:,:,iz,ii);
			Dflt(:,:,iz,ii) = bsxfun(@times, d, flt);               % filtered k-space
		end
	end
end

Dflt = reshape(Dflt, sz);

return;


function Dflt = sub_test()

zshape = [zeros(1,6) ones(1,12) zeros(1,6)];   % 24 slices
N = 128;
for ii = 1:length(zshape)
	D(:,:,ii) = zshape(ii)*phantom(N);
end

fltRadius = N/10;
transitionWidth = fltRadius/4;
Dflt = toppe.utils.imfltfermi(D, fltRadius, transitionWidth, 'rect');

return;
