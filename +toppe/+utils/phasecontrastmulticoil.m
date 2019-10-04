function pc = phasecontrastmulticoil(ims1, ims2)
% function pc = phasecontrastmulticoil(ims1, ims2)
%
% Calculate phase contrast between two multi-coil 2D/3D images.
%
% size(ims1) = nx x ny x ncoils OR nx x ny x nz x ncoils
%
% OTPUT:
%   pc      
%

if ndims(ims1) > 3
	[nx, ny, nz, ncoils] = size(ims1);
else
	[nx, ny, ncoils] = size(ims1);
	nz = 1;
end

pc = zeros(nx, ny, nz);

for coil = 1 : ncoils
	if ndims(ims1) > 3
		im1 = ims1(:,:,:,coil);
		im2 = ims2(:,:,:,coil);
	else
		im1 = ims1(:,:,coil);
		im2 = ims2(:,:,coil);
	end
	pc = pc + im1.*conj(im2);
end;

pc = angle(pc);

return;
