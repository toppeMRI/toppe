function pc = phasecontrast3d(ims1, ims2)
% Calculate phase contrast between two multi-coil 3D images.
%
% function pc = phasecontrast3d(ims1, ims2)
%
% size(ims1) = nx x ny x nz x ncoils
%
% OTPUT:
%   pc      
%
% $Id: phasecontrast3d.m,v 1.1 2018/11/02 14:38:31 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/phasecontrast3d.m,v $

[nx, ny, nz, ncoils] = size(ims1);

pc = zeros(nx, ny, nz);

for coil = 1 : ncoils
	im1 = ims1(:,:,:,coil);
	im2 = ims2(:,:,:,coil);
	pc = pc + im1.*conj(im2);
end;

pc = angle(pc);

return;
