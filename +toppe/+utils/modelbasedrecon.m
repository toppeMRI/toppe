function [ims] = modelbasedrecon(kspace, data, imsize, fov, varargin)
% Performs a model-based iterative reconstruction using a quadratic
% penalized weighted least squares objective function with preconditioned
% conjugate gradient, specifically using MIRT's qpwls_pcg1 function.
%
% Needs MIRT.
%
% Can recon multiple images with the same A (such as a 2D acquisition over
% slices, or fMRI timepoints with the same 3D system matrix).
%
% function [im] = modelbasedrecon(A, data, varargin)
%
% Inputs:
%   kspace      [ndat nshot 2|3]   input kspace
%   data        [ndat nshot ncoil nim]   k-space data
%   imsize      [2 1] or [3 1]     size of output image
%   fov         [2 1] or [3 1]     fov of output image

% Outputs:
%   ims         [nx ny | nz]


% Optional inputs:
%   sensemaps   [nx ny nz ncoil]
%   other varargin variables - see defaults below

% Melissa Haskell
% July 2020

if nargin == 1 && streq(kspace, 'test') 
    modelbasedrecon_test;    % 2d 
    modelbasedrecon_test3d;  % 3d
    return
end


%% Hardcoded TOPPE variables
dt = 4e-6; % toppe sampling time is 4 us.


%% make sure MIRT is installed
% (note: not sure this is the best way to check this)
if exist('qpwls_pcg1','file') ~= 2
    error(['Error: modelbasedrecon.m requires the Michigan Image ',...
        'Reconstruction Toolbox to run qpwls_pcg1. '])
end


%% check inputs

% compare dimensions of kspace & input data
ndat = size(kspace,1);
nshot = size(kspace,2);
if numel(size(data)) < 3
    ncoil = 1; nim = 1;
elseif numel(size(data)) < 4
    ncoil = size(data,3); nim = 1;
else
    ncoil = size(data,3); nim = size(data,4);
end
if ndat ~= size(data,1), error('Incorrect variable input sizes.'); end
if nshot ~= size(data,2), error('Incorrect variable input sizes.'); end

% make sure image space info same dimentions
if size(imsize) ~= size(fov), error('Incorrect variable input sizes.'); end

%% Set defaults and parse varargin

% set defaults
arg.zmap          = [];
arg.useParallel   = false;
arg.np            = round(.75*(feature('numcores')));  % number of parallel workers
arg.quiet         = true;
arg.W             = 1;  % default data weights
arg.C             = 0;  % default quadratic penalty (i.e. dafult is no reg)
arg.x0            = zeros(prod(imsize),nim); % initial image estimate
arg.niter         = 5;  % default number of iteration of pcg
arg.stop_diff_tol = 1e-3; % default stopping tolerance
arg.sensemaps     = ones(imsize); % if none given, assume one uniform coil

% parse variable input arguments
arg = toppe.utils.vararg_pair(arg, varargin);

% NOTE! qpwls_pcg1 has many other options, but to use those their defaults
% would need to be added to this section of modelbasedrecon.m

% nufft options and L (taken from TOPPE's reconSoS function
nx = imsize(1); ny = imsize(2);
nufft_args = {[nx,ny],[6,6],[2*nx,2*ny],[nx/2,ny/2],'minmax:kb'};
% mask = true(Ny,Nx); % Mask for support
L = 6;

%% Setup image geometry & timing, then create forward model fatrix

% reshape kspace and raw data by concatenating shots
kspace_cat = reshape(kspace, [ndat*nshot, size(kspace,3)]);
data = reshape(data,[ndat*nshot*ncoil, nim]);

% create timing vector
t = repmat(0:dt:dt*(ndat-1), [1 nshot]);

% setup either 2d or 3d fatrix object for image recon
if numel(imsize) == 2
    mirt_image_geom = image_geom('nx', nx, 'ny', ny, 'fov', fov);
    if ~isempty(arg.zmap)
        mirt_output = evalc(strcat("A0=Gmri(kspace_cat(:,1:2), ",...
            "mirt_image_geom.mask,'fov', mirt_image_geom.fov, 'ti', t, ",...
            "'zmap', 1i*2*pi*arg.zmap, 'L', L,'nufft',nufft_args)"));
    else
        mirt_output = evalc(strcat("A0=Gmri(kspace_cat(:,1:2), ",...
            "mirt_image_geom.mask,'fov', mirt_image_geom.fov, 'ti', t)"));
    end
    A = Asense(A0, arg.sensemaps);
else
    %%%%% 7/13/20 note - 3d hasn't been tested yet
    nz = imsize(3);
    nufft_args = {[nx,ny,nz],[6,6,6],[2*nx,2*ny,2*nz],[nx/2,ny/2,nz/2],'minmax:kb'};
    mask = true(nx,ny,nz); % Mask for support
    L = 6;
    mirt_image_geom = image_geom('nx', nx, 'ny', ny, 'nz', nz, ...
        'fov', fov);
    if ~isempty(arg.zmap)
        mirt_output = evalc(strcat("A0 = Gmri(kspace_cat, mirt_image_geom.mask,",...
            "'fov', mirt_image_geom.fov, 'ti', t, 'zmap', 1i*2*pi*arg.zmap,", ...
            "'L', L, 'nufft', nufft_args)"));
    else
        mirt_output = evalc(strcat("A0 = Gmri(kspace_cat, mirt_image_geom.mask,",...
            "'fov', fov, ", ...
            "'L', L, 'nufft', nufft_args)"));
    end
    A = Asense(A0, arg.sensemaps);
end


%% call qpwls_pcg1 from MIRT

% initialize 2d matrix to hold image vectors
npix = prod(imsize);
ims_vec_all = zeros(npix,nim);

% get image recon input parameters
x0 = arg.x0; W = arg.W; C = arg.C; niter = arg.niter;
stop_diff_tol = arg.stop_diff_tol;

% call either in parallel or serially
if arg.useParallel
    if ~arg.quiet; fprintf('Reconstructing in parallel... \n'); end
    
    % start parallel pool if needed
    p = gcp('nocreate');
    if isempty(p)
        parpool(arg.np)
    end   
    
    % loop over images and reconstruct
    parfor ii = 1:nim
        [ims_vec_all(:,ii), ~] = ...
            qpwls_pcg1(x0(:,ii), A, W, data(:,ii), C, ...
            'niter', niter, 'stop_diff_tol', stop_diff_tol,...
            'isave', 'last');
    end
    
    if ~arg.quiet; fprintf('done.\n'); end
else
    if ~arg.quiet; fprintf('Reconstructing... \n'); end
    
    % loop over images and reconstruct
    for ii = 1:nim
        mirt_output = evalc(strcat("[ims_vec_all(:,ii), ~] = ",...
            "qpwls_pcg1(arg.x0(:,ii), A, W, data(:,ii), C, ",...
            "'niter', niter, 'stop_diff_tol', stop_diff_tol,",...
            "'isave', 'last');"));
        
    end
    
    if ~arg.quiet; fprintf('done.\n'); end
end

%% Return all images in proper sized matrices
if numel(imsize) == 2
    ims = reshape(ims_vec_all, [nx ny nim]);
else
    ims = reshape(ims_vec_all, [nx ny nz nim]);
end


end




function [] = modelbasedrecon_test
% basic 2D test 

import toppe.*
import toppe.utils.*
imfig = 10;

% image parameters
n = 128;        % image size 
fov = 200;      % field of view (mm)
pe_dir = 1;     % phase encode direction
pad = 10;
ncoil = 4;

%% 2d test

% create image geometry & im, display
ig = image_geom('nx', n, 'fov', fov); 
xtrue = padarray(phantom('Modified Shepp-Logan',n-2*pad),[pad pad], 0);
figure(imfig); subplot(131); imagesc(abs(xtrue)); title('ground truth')
axis image; colormap gray; colorbar

% Synthesize k-space data 
f.traj = 'cartesian';
[kspacefull, ~, ~] = mri_trajectory(f.traj, {}, [n n], ig.fov);
smaps = ir_mri_sensemap_sim('nx',n,'ny',n,'ncoil',ncoil);
A0 = Gmri(kspacefull, ig.mask,'fov', ig.fov);
A = Asense(A0, smaps);
dattestfull = A*xtrue(:);

% Undersample
kmask = true(n,n);
for iy = 1:n
    kmask(round((mod(iy,2)+1):2:end),iy) = false; % R=2 CAIPI
end
%r = ((-n/8):(n/8)) + n/2;
%kmask(r,r) = true; % densely sampled center
dattestfull = reshape(dattestfull, [n n ncoil]);
kspacefull = reshape(kspacefull, [n n, 2]);
for ic = 1:ncoil
    dtmp = dattestfull(:,:,ic);
    dattest(:, ic) = dtmp(kmask);
end
for id = 1:2
    ktmp = kspacefull(:, :, id);
    kspace(:, id) = ktmp(kmask);
end

% recon and display
kspace = reshape(kspace,[size(kspace,1),1,2]);
dattest = reshape(dattest,[size(dattest,1),1,ncoil]);
[ims] = modelbasedrecon(kspace, dattest, [n n], [fov fov], ...
    'sensemaps', smaps, ...
    'niter', 20);
figure(imfig); subplot(132); imagesc(abs(ims)); title('MBIR')
axis image; colormap gray; colorbar
disp('2d single slice test complete.')

%subplot(133); im(kmask);
%colormap default

%% stack of 2d test
nim=5;
dattest = repmat(dattest,[1 1 1 nim]);
[ims] = modelbasedrecon(kspace, dattest, [n n], [fov fov], ...
    'sensemaps', smaps, ...
    'niter', 20);
figure(imfig); subplot(132); imagesc(abs(ims(:,:,1))); title('MBIR 1st im')
axis image; colormap gray; colorbar
subplot(133); imagesc(abs(ims(:,:,end))); title('MBIR last im')
axis image; colormap gray; colorbar
disp('2d multislice test complete.')

%% stack of 2d test parallel
try
nim=5;
dattest = repmat(dattest,[1 1 1 nim]);
[ims] = modelbasedrecon(kspace, dattest, [n n], [fov fov], 'sensemaps', ...
    smaps, 'useParallel', true);
%figure(imfig); subplot(132); imagesc(abs(ims(:,:,1))); title('MBIR 1st im')
axis image; colormap gray; colorbar
subplot(133); imagesc(abs(ims(:,:,end))); title('MBIR last im')
axis image; colormap gray; colorbar
disp('2d multislice parallel compute test complete.')
catch ME
    warning('Parallel test failed');
end

end




function [] = modelbasedrecon_test3d
% basic 3D test 

fprintf('3d test...');

import toppe.*
import toppe.utils.*

imfig = 11;

nd = 3;   % 3 dimensions

% image parameters
N = [64 64 16];        % image size 
[nx ny nz] = deal(N(1), N(2), N(3));
dx = 200/nx;   % voxel size (mm)
FOV = dx*N;      % field of view (mm)
pad = 10;
ncoil = 4;

% object
%xtmp = padarray(phantom('Modified Shepp-Logan',nx-2*pad),[pad pad], 0);
xtmp = phantom('Modified Shepp-Logan',nx);
xtrue = zeros(N);
xtrue(:,:,3) = xtmp;
for iz = 4:(nz-2)
    xtrue(:,:,iz) = xtrue(:,:,iz-1)';
end
figure(imfig); subplot(131); im(abs(xtrue(:,:,end/2))); title('ground truth (middle slice)')

% Synthesize k-space data 
f.traj = 'cartesian';
[kspacefull, ~, ~] = mri_trajectory(f.traj, {}, N, FOV);  % [prod(N) 3]
smaps = ir_mri_sensemap_sim('nx', nx, 'ny', ny, 'nz', nz, 'ncoil', ncoil);  % [[N] ncoil]
nufft_args = {N, [6,6,6], [2*nx, 2*ny, 2*nz], [nx/2, ny/2, nz/2], 'minmax:kb'};
A0 = Gmri(kspacefull, true(N), 'fov', FOV, ...
    'nufft', nufft_args);
A = Asense(A0, smaps);
dattestfull = A*xtrue(:);

% Undersample
kmask = false(N);
for iy = 1:ny
    kmask(round((mod(iy,2)+1):2:end),iy,:) = true; % R=2 CAIPI
end
%rx = ((-nx/8):(nx/8)) + nx/2;
%rz = ((-nz/8):(nz/8)) + nz/2;
%kmask(rx,rx,nz) = true; % densely sampled center
%R = sum(kmask(:))/prod(N)  % undersampling factor
dattestfull = reshape(dattestfull, [nx ny nz ncoil]);
kspacefull = reshape(kspacefull, [nx ny nz nd]);
for ic = 1:ncoil
    dtmp = dattestfull(:,:,:,ic);
    dattest(:, ic) = dtmp(kmask);
end
for id = 1:nd
    ktmp = kspacefull(:, :, :, id);
    kspace(:, id) = ktmp(kmask);
end

% recon and display
kspace = reshape(kspace,[size(kspace,1), 1, nd]);
dattest = reshape(dattest,[size(dattest,1),1,ncoil]);
[ims] = modelbasedrecon(kspace, dattest, N, FOV, ...
    'sensemaps', smaps, ...
    'niter', 10);
figure(imfig); subplot(132); im(abs(ims(:,:,end/2))); title('MBIR (middle slice)')
figure(imfig); subplot(133); im(abs(ims)); title('MBIR')
axis image; colormap default; colorbar

fprintf('\t3d test copmlete\n');

end

