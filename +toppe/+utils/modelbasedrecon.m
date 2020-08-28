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
        A0 = Gmri(kspace_cat(:,1:2), mirt_image_geom.mask,'fov',...
            mirt_image_geom.fov, 'ti', t, 'zmap', 1i*2*pi*arg.zmap,...
            'L', L,'nufft',nufft_args);
    else
        A0 = Gmri(kspace_cat(:,1:2), mirt_image_geom.mask,'fov',...
            mirt_image_geom.fov, 'ti', t);
    end
    A = Asense(A0, arg.sensemaps);

else
    %%%%% 7/13/20 note - 3d hasn't been tested yet
    nz = imasize(3);
    mirt_image_geom = image_geom('nx', nx, 'ny', ny, 'nz', nz, ...
        'fov', fov);
    if ~isempty(arg.zmap)
        A0 = Gmri(kspace_cat, mirt_image_geom.mask,'fov',...
            mirt_image_geom.fov, 'ti', t, 'zmap', 1i*2*pi*arg.zmap, ...
            'L', L, 'nufft', nufft_args);
    else
        A0 = Gmri(kspace_cat, mirt_image_geom.mask,'fov',...
            mirt_image_geom.fov, 'ti', t);
    end
    A = Asense(A0, arg.sensemaps);
end


%% call qpwls_pcg1 from MIRT

% initialize 2d matrix to hold image vectors 
npix = prod(imsize);
ims_vec_all = zeros(npix,nim);

% call either in parallel or serially
if arg.useParallel
    if ~arg.quiet; fprintf('Reconstructing in parallel... \n'); end
   
    % start parallel pool if needed
    p = gcp('nocreate');
    if isempty(p)
        parpool(arg.np)
    end
    
    % get image recon input parameters
    x0 = arg.x0; W = arg.W; C = arg.C; niter = arg.niter;
    stop_diff_tol = arg.stop_diff_tol;
    
    % loop over images and reconstruct
    parfor ii = 1:nim
        [ims_vec_all(:,ii), ~] = qpwls_pcg1(x0(:,ii), A, W,...
            data(:,ii), C, 'niter', niter, 'stop_diff_tol', stop_diff_tol);
    end
    
    if ~arg.quiet; fprintf('done.\n'); end
else
    if ~arg.quiet; fprintf('Reconstructing... \n'); end
    
    % loop over images and reconstruct
    for ii = 1:nim
        [ims_vec_all(:,ii), ~] = qpwls_pcg1(arg.x0(:,ii), A, arg.W,...
            data(:,ii), arg.C, 'niter', 10, 'stop_diff_tol', 1e-3);
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
