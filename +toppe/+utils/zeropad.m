% zero fill 3D images, (optionally) Fermi-filtered (in-plane) to reduce ringing
%
% function imout = zeropad(imin, res, transitionwidth, fltgeom)
%
% imin               [nx ny nz]
% res                [nxtarget nytarget nztarget]      
% transitionwidth    Fermi filter transition width (pixels). Default: 2.
% fltgeom               'circ' (default) or 'rect'

function imout = zeropad(imin, res, transitionwidth, fltgeom)

import toppe.utils.*

if nargin < 3
	dofilter = false;
else
	dofilter = true;
end

if ~exist('fltgeom', 'var')
	fltgeom = 'circ';
end

[nx ny nz] = size(imin);
if nx ~= ny
	error('image must be square in x/y');
end

% zero-pad
dout = zeros(res);
rangex = (res(1)/2+1-nx/2):(res(1)/2+1+nx/2-1);
rangez = (res(3)/2+1-nz/2):(res(3)/2+1+nz/2-1);
d = cfftn(imin, 'forward');
dout(rangex, rangex, rangez) = d;

% filter
if dofilter
	if ~exist('transitionwidth', 'var')
		transitionwidth = 2;    % filter transition width (pixels)
	end
	flt = fermi2d(res(1), size(d,1), transitionwidth, fltgeom);
	for iz = 1:size(dout,3)
		dout(:,:,iz) = bsxfun(@times, dout(:,:,iz), flt);
	end
end

imout = cfftn(dout, 'inverse');

return

