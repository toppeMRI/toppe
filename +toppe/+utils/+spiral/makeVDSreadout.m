function [gx_tmp, gy_tmp, gz_tmp, paramsint16, paramsfloat] = ...
    makeVDSreadout(fov, F0, F1, F2, npix, nshots, maxs_xy, ...
    spiraldir, system, varargin)
% Wrapper for 2D balanced VDS encoding. Optional multishot.

% INPUTS:
%  fov - imaging field of view, cm
%  F0 - initial field of view for k-space acquisition, cm
%  F1 - F1 linear coefficient for vds
%  F2 - F2 quadratic coefficient for vds
%  npix   - number of pixels for matrix
%  nshots - number of shots
%  maxs_xy - max x, y slew in unit of G/cm/ms
%  spiraldir - spiral direction (out=1, in=2, inout = 3)

% Options:
%   crushers    [ncrush 3] Gradient crushers to add to end of readout
%   spfrac      [1 1]      Spiral fraction - ratio of spiral in to spiral 
%                           out for spiral in - out sequence
%   nFIDstart   [1 1]      # of times to sample k=0 at beginning of readout
%   nFIDend     [1 1]      # of times to sample k=0 at end of readout
%   gmax        [1 1]      peak gradient (G/cm). Default: 5 G/cm

% inital version by Amos Cao
% updated June 2020 by Melissa Haskell, mhask@umich.edu


%% parse inputs and set defaults

% set  defaults
arg.spfrac    = 0.2;
arg.nFIDstart = 0; 
arg.nFIDend   = 0;  
arg.crushers  = [];

% system defaults
arg.gmax = 5;	 % G/cm

% parse inputs
arg = vararg_pair(arg, varargin);

% toppe sampling time
T = 4e-6;	 % Seconds

% balancing optimization defaults
maxiter_bal_in = 50;
maxiter_bal_out = 50;


%% Generate VDS for gx and gy
smax = maxs_xy*1000; % multiply by 1000 to get to T/m/s

rmax = npix/(2*fov);		% cm^(-1)
N = nshots;

disp('Calculating gradients...');

%%%% Main spiral out - will flip if spiral in indicated
[~,gout,~,~,~,~] = toppe.utils.spiral.mintgrad.vds(smax,arg.gmax,T,nshots,F0,F1,F2,rmax);
nsamp2 = length(gout);
gx_out = real(gout(:));
gy_out = imag(gout(:));

% mini spiral out that will become spiral in for spiral in-out
if spiraldir == 3
    [~,gin,~,~,~,~] = toppe.utils.spiral.mintgrad.vds(smax,arg.gmax,T,nshots,F0,F1,F2,rmax*arg.spfrac);
    nsamp1 = length(gin);
    gx_in = real(gin(:));
    gy_in = imag(gin(:));
else
    nsamp1 = 0;
end




%% Balance spirals
fprintf('Creating balancers...');

% Balance the spiral out, calling toppe.utils.gbalance multiple times if
% needed. 
iter = 1;
gx_tmp = gx_out; gy_tmp = gy_out;
while max(abs([trapz(gx_tmp) trapz(gy_tmp)])) > 1e-6 && iter < maxiter_bal_out
    
    gxb = toppe.utils.gbalance(gx_tmp,maxs_xy,0); fprintf('.');
    gyb = toppe.utils.gbalance(gy_tmp,maxs_xy,0); fprintf('.');
    lengthdiff = length(gxb)-length(gyb);
    
    % Rebalance the shorter waveform with the extra time we have
    if lengthdiff > 0
        gyb = toppe.utils.gbalance(gy_tmp,maxs_xy,length(gxb));
    elseif lengthdiff < 0
        gxb = toppe.utils.gbalance(gx_tmp,maxs_xy,length(gyb));
    end
    
    gxbal = [gx_tmp; gxb];
    gybal = [gy_tmp; gyb];
    
    gx_tmp = gxbal;
    gy_tmp = gybal;
    
    iter = iter + 1;
end
gx_out_bal = gx_tmp; gy_out_bal = gy_tmp;
nBalance_out = numel(gx_out_bal) - numel(gx_out(:));

% balance spiral in if we have a spiral in-out sequence
if spiraldir == 3
    gx_tmp = gx_in; gy_tmp = gy_in; iter = 1;
    while max(abs([trapz(gx_tmp) trapz(gy_tmp)])) > 1e-6 && iter < maxiter_bal_in
        
        gxb = toppe.utils.gbalance(gx_tmp,maxs_xy,0); fprintf('.');
        gyb = toppe.utils.gbalance(gy_tmp,maxs_xy,0); fprintf('.');
        lengthdiff = length(gxb)-length(gyb);
        
        % Rebalance the shorter waveform with the extra time we have
        if lengthdiff > 0
            gyb = toppe.utils.gbalance(gy_tmp,maxs_xy,length(gxb));
        elseif lengthdiff < 0
            gxb = toppe.utils.gbalance(gx_tmp,maxs_xy,length(gyb));
        end
        
        gxbal = [gx_tmp; gxb];
        gybal = [gy_tmp; gyb];
        
        gx_tmp = gxbal;
        gy_tmp = gybal;
        
        iter = iter + 1;
    end
    gx_mini_in_bal = gx_tmp; gy_mini_in_bal = gy_tmp;
    nBalance_in = numel(gx_mini_in_bal) - numel(gx_in);
else
    gx_mini_in_bal = []; gy_mini_in_bal = [];
    nBalance_in = 0;
end


% Check to make sure we're actually balanced
if max(abs([trapz(gx_out_bal) trapz(gy_out_bal)])) > 1e-6
    fprintf('K ending: %0.7f,%0.7f\n',trapz(gx_out_bal),trapz(gy_out_bal));
    warn('Failed to balance sprial out, proceeding...');
end
if max(abs([trapz(gx_mini_in_bal) trapz(gy_mini_in_bal)])) > 1e-6
    fprintf('K ending: %0.7f,%0.7f\n',trapz(gx_mini_in_bal),trapz(gy_mini_in_bal));
    warn('Failed to balance sprial in, proceeding...');
end

fprintf(' done.\n');


%% Combine spiral in and out gradients 
gx_tmp = [gx_mini_in_bal(end:-1:1); gx_out_bal];
gy_tmp = [gy_mini_in_bal(end:-1:1); gy_out_bal];
gz_tmp = zeros(size(gx_tmp));


%% Additional sequence modifications and checks

% Do multishot rotations
if nshots > 1
    [gx_tmp, gy_tmp] = shotrot(gx_tmp,gy_tmp,nshots);
    gz_tmp = zeros(size(gx_tmp));
end

% Flip if spiral in (can simply flip bc balanced)
if spiraldir == 2
    gx_tmp = flipud(gx_tmp);
    gy_tmp = flipud(gy_tmp);
    gz_tmp = flipud(gz_tmp);
    
    nBalance_in = nBalance_out;
    nBalance_out = 0;
end

% Add FID points at the beginning & end (e.g. to calculate bulk B0 shift per TR)
gx_tmp = [zeros(arg.nFIDstart,size(gx_tmp,2)); gx_tmp; zeros(arg.nFIDend,size(gx_tmp,2))];
gy_tmp = [zeros(arg.nFIDstart,size(gy_tmp,2)); gy_tmp; zeros(arg.nFIDend,size(gy_tmp,2))];
gz_tmp = [zeros(arg.nFIDstart,size(gz_tmp,2)); gz_tmp; zeros(arg.nFIDend,size(gz_tmp,2))];

% check for non-zeros entries at first and last indices
if sum(gx_tmp(1,:)~=0) || sum(gy_tmp(1,:) ~=0)
    gx_tmp = [zeros(4,size(gx_tmp,2)); gx_tmp];
    gy_tmp = [zeros(4,size(gy_tmp,2)); gy_tmp];
    gz_tmp = [zeros(4,size(gz_tmp,2)); gz_tmp];
    arg.nFIDstart = 4;
end
if sum(gx_tmp(end,:)~=0) || sum(gy_tmp(end,:) ~=0)
    gx_tmp = [gx_tmp; zeros(4,size(gx_tmp,2))];
    gy_tmp = [gy_tmp; zeros(4,size(gy_tmp,2))];
    gz_tmp = [gz_tmp; zeros(4,size(gz_tmp,2))];
    arg.nFIDend = 4;
end

gx = gx_tmp;
gy = gy_tmp;
gz = gz_tmp;

% add crushers if specified in variable inputs
if ~isempty(arg.crushers)
    if numel(size(arg.crushers)) == 3  % needed for multishot
        gx = [gx; arg.crushers(:,:,1)];
        gy = [gy; arg.crushers(:,:,2)];
        gz = [gz; arg.crushers(:,:,3)];
    else
        gx = [gx; arg.crushers(:,1)];
        gy = [gy; arg.crushers(:,2)];
        gz = [gz; arg.crushers(:,3)];
    end
    nCrusher = size(arg.crushers,1);
else
    nCrusher = 0;
end

% size check
if any(size(gx) ~= size(gy)) || any(size(gx) ~= size(gz)) || any(size(gz) ~= size(gy))
    error('Gx gy and gz aren''t the same length... this should not happen!');
end


%% Store parameters and write readout.mod file

% Store parameters that are integers
paramsint16(1) = npix;
paramsint16(2) = nsamp2;
paramsint16(3) = nshots;
paramsint16(6) = arg.nFIDstart;
paramsint16(7) = arg.nFIDend;
paramsint16(8) = spiraldir;
paramsint16(9) = nBalance_out;
paramsint16(10) = nsamp1;
paramsint16(11) = nBalance_in;
paramsint16(12) = fov; 
paramsint16(13) = nCrusher;

% Store parameters that are floats
paramsfloat(1) = smax;
paramsfloat(2) = arg.gmax;
paramsfloat(3) = T;
paramsfloat(4) = N;
paramsfloat(5) = F0;
paramsfloat(6) = F1;
paramsfloat(7) = F2;
paramsfloat(8) = rmax;
paramsfloat(9) = fov; 
paramsfloat(10) = F1;
paramsfloat(11) = F2;

toppe.writemod(system,'gx',gx,'gy',gy,'gz',gz_tmp,'ofname','readout.mod','hdrints',paramsint16,'hdrfloats',paramsfloat);

end


function [gxrot, gyrot] = shotrot(gx,gy,nshots)
%shotrot - Rotates gradient shots evenly for the specified number of shots
% Inputs should be nsamp x nwaveforms complex vectors for x and y gradients or kspace

% Author: Amos Cao, M.S.E., Biomedical Engineering
% Functional MRI Laboratory, University of Michigan
% Email address: amoscao@umich.edu
% Feb 2019

shotrot = 2*pi/nshots;
rotvec = exp(-1i*(0:nshots-1)*shotrot);
rotvec = repmat(rotvec,size(gx,2),1); rotvec = rotvec(:).'; %Repeat each entry x number of waveforms inputted
g = gx + 1i*gy;
gmulti = repmat(g,1,nshots);
rotmatrix = repmat(rotvec,length(g),1);
grot = gmulti.*rotmatrix; % Rotate
gxrot = real(grot);
gyrot = imag(grot);
end

