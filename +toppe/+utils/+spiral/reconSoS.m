function [imsos ims dcf] = reconSoS(dat, kx, ky, fov, imsize, varargin)
% Reconstruct fully-sampled stack-of-spirals/stars data via adjoint nufft.
% With (optional) fieldmap correction. Needs MIRT.
%
% function [ims imsos] = reconSoS(dat, kx, ky, fov, imsize, varargin)
%
% Inputs:
%   dat         [ndat nleafs nz ntp ncoils]        acquired (kspace) data [ndat nleafs nz ntp ncoils]
%   kx          [ndat leafs]                       readout trajectory (cycles/cm)
%                                                 No support for timecourse
%                                                 varying readout yet (ie
%                                                 rotating spiral)
%   ky          [ndat leafs]                       ""
%   fov         [2 1] or [3 1]                     FOV in x/y/(z) (cm). z FOV needed if fieldmap is provided.
%   imsize      [nx ny]                            reconstructed image matrix size
% Options:
%   dcf         density compensation function
%   zmap        [imsize nz] 3D image volume with field map values (Hz). Matrix size must match imsize. See Gmri.m.
%   ti          sample times
%   useParallel do recon using parpool (requires parallel toolbox)
%
% $Id: reconSoS.m,v 1.8 2018/11/12 13:38:46 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+spiral/reconSoS.m,v $

import toppe.utils.spiral.*

% Set defaults and parse varargin
arg.zmap         = [];
arg.ti           = [];
arg.spatialShift = [0 0];
arg.dcf          = [];
arg.useParallel = false;
arg = vararg_pair(arg, varargin);

[ndat nleafs nz ntp ncoils] = size(dat);

% calculate density compensation function (if not provided)
if isempty(arg.dcf)
    knorm = bsxfun(@rdivide, [fov(1)*kx(:) fov(2)*ky(:)], imsize(1:2)); % normalize to [-0.5 0.5]
    dcf = mdcf_voronoi(knorm);
else
    dcf = arg.dcf;
end

% system matrix options
npix = imsize(1);
%xinit = zeros(npix^2,1);
Ny = npix;
Nx = npix;
nufft_args = {[Ny,Nx],[6,6],[2*Ny,2*Nx],[Ny/2,Nx/2],'table',2^11,'minmax:kb'};
mask = logical(ones(Ny,Nx)); % Mask for support
L = 6;

% Initialize system matrix for each z location
for iz = 1:nz
    if ~isempty(arg.zmap)    % do nufft with fieldmap correction
        if isempty(arg.ti)
            error('You must provide ti (sample times) as well as zmap');
        end
        
        % Initialize a Gmri object for every slice location
        A{iz} = Gmri([fov(1)*kx(:) fov(2)*ky(:)], mask, ...
            'ti', arg.ti, 'zmap', 1i*2*pi*arg.zmap(:,:,iz), 'L', L, 'nufft', nufft_args);
    else
        A{iz} = Gmri([fov(1)*kx(:) fov(2)*ky(:)],mask,'nufft',nufft_args);
    end
end

% Do IFT along kz, then in-plane nufft per slice and per timepoint
dat = sub_iftz(dat);
ims = zeros([imsize nz ntp ncoils]);
if arg.useParallel
    fprintf('Reconstructing in parallel... ');
    parfor icoil = 1:ncoils
        for itp = 1:ntp
            for iz = 1:nz
                d2d = squeeze(dat(:,:,iz,itp,icoil));   % [ndat nleafs]
                imtmp = A{iz}'*(d2d(:).*dcf);
                ims(:,:,iz,itp,icoil) = reshape(imtmp, imsize(1:2));
            end
        end
    end
    fprintf('done.\n');
else
    toppe.utils.textprogressbar('Reconstructing: ');
    for icoil = 1:ncoils
        for itp = 1:ntp
            toppe.utils.textprogressbar(100*(itp+((icoil-1)*ntp))/(ntp*ncoils));
            for iz = 1:nz
                d2d = squeeze(dat(:,:,iz,itp,icoil));   % [ndat nleafs]
                imtmp = A{iz}'*(d2d(:).*dcf);
                ims(:,:,iz,itp,icoil) = reshape(imtmp, imsize(1:2));
            end
        end
    end
    toppe.utils.textprogressbar(' done.');
end
imsos = sqrt(sum(abs(ims).^2,5)); % Do root-sum-of-squares coil combination

return;

function x = sub_iftz(y)
n = size(y,3);
n = 1;
x = fftshift(ifft(ifftshift(y,3), [], 3), 3 )*n;
return;
