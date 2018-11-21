function dat = preprocessSoS(dat, kx, ky, varargin)
% Apply data shift, spatial shift, b0 eddy current demodulation, to stack-of-spirals/stars data.
%
% Inputs:
%  dat      [ndat nleafs nz ncoils]             Acquired (raw) kspace data
%  kx       [ndat leafs] (cycles/cm)            kx space
%  ky       [ndat leafs] (cycles/cm)            ky space
% Options:
%  datShift       [1]        Acquisition/gradient delay (samples)
%  spatialShift   [dx dy]    In-plane shift (mm)
%  b0eddy         [ndat 1]   Measured B0 eddy current.
%
% $Id: preprocessSoS.m,v 1.4 2018/11/11 16:29:36 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/+utils/+spiral/preprocessSoS.m,v $

% Set defaults and parse varargin
arg.datShift     = 0.0;                      % shift data frame by this many samples (to correct for gradient/acquisition delay)
arg.spatialShift = [0 0];                    % in-plane (x, y) spatial shift (mm)
arg.b0eddy       = [];
arg.offres       = 0;                        % Global off-resonance offset (Hz)
arg.flipz        = false;
arg = vararg_pair(arg, varargin);
 
% shift data frame (grad/acq delay)
[ndattmp nl nz ncoils] = size(dat);
[a b c] = ndgrid( (1:ndattmp)+arg.datShift, 1:nl, 1:nz);
for ic = 1:ncoils
	dat(:,:,:,ic) = interpn(dat(:,:,:,ic), a, b, c);
end

% flip in z
if arg.flipz
	dat = flipdim(dat, 3);
end

[ndat nl nz ncoils] = size(dat);

% b0 eddy demodulation
if ~isempty(arg.b0eddy)
	for ic = 1:ncoils
		for iz = 1:nz
			dat(:,:,iz,ic) = exp(-1i*arg.b0eddy/180*pi).*dat(:,:,iz,ic);
		end
	end
end

% apply spatial shift (in-plane)
[dx dy] = deal(arg.spatialShift(1)/10, arg.spatialShift(2)/10);  % cm
for ic = 1:ncoils
	for iz = 1:nz
		dat(:,:,iz,ic) = dat(:,:,iz,ic).*exp(1i*2*pi*kx*dx).*exp(1i*2*pi*ky*dy);
	end
end

