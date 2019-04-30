function [ims imsos d]= epi_recon(pfile, readoutfile, gradDelay, phsOffset)
% Reconstruct 2D EPI data acquired with ISMRM2019 "live" demo
%
% Inputs:
%  pfile
%  readoutfile   Default: 'readout.mod'
%  gradDelay     gradient/acquisition delay (sec)
%
% Output:
%  ims           [nx ny ncoils]    
%  imsos         coil-combined (root-sum-of-squares) image
%  d             raw (k-space) data

%addpath ~/gitlab/toppe/
%import toppe.*
%import toppe.utils.*

if ~exist('readoutfile','var')
	readoutfile = 'readout.mod';
end
if ~exist('gradDelay','var')
	gradDelay = 0;
end
if ~exist('phsOffset','var')
	phsOffset = 0;
end

%% get readout file header
[~,gx,~,~,~,hdrints] = toppe.readmod(readoutfile);
ndat = size(gx,1);
N    = hdrints(1);       % image size
nes  = hdrints(3);       % echo spacing (number of 4us samples)
npre = hdrints(4);       % number of samples before start of readout plateau of first echo
nshots = size(gx,2);

%% load raw data
d = toppe.utils.loadpfile(pfile); %, 1, 2, 2);               % int16, size [ndat ncoils nslices nechoes nviews] = [ndat ncoils 1 1 nshots]
d = permute(d,[1 5 2 3 4]);         % [ndat nshots ncoils].
d = double(d);
d = flipdim(d,1);        % data is stored in reverse order (for some reason)
[ndat nshots ncoils] = size(d);

%% apply gradient/acquisition delay and phase offset
d = circshift(d, 1);
d(:,2:2:end,:) = exp(-1i*phsOffset)*d(:,2:2:end,:);
nup = 10;       % upsampling factor
dt = 4e-6;      % raster time (sec)
nshift = round(gradDelay/dt*nup);
dup = interpft(d, nup*ndat, 1);
dup = circshift(dup,nshift,1);
%d = dup(1:nup:end,:,:);
for ii = 1:nshots
	for jj = 1:ncoils
		d(:,ii,jj) = decimate(dup(:,ii,jj),nup);
	end
end

%% sort data into 2D NxN matrix
d2d = zeros(N,N,ncoils);
etl = N/nshots;    % echo-train length
for ic = 1:ncoils
	for ii = 1:nshots
		for jj = 1:etl
			istart = npre + (jj-1)*nes + 1;
			dtmp = d(istart:(istart+N-1), ii, ic);  % one echo

			% flip echo as needed 
			if mod(ii,2)
				if mod(jj-1,2)
					dtmp = flipdim(dtmp,1);           
				end
			else
				if mod(jj,2)
					dtmp = flipdim(dtmp,1);
				end
			end

			% phase-encode index
			iy = (jj-1)*nshots + ii;

			d2d(:,iy,ic) = dtmp;
		end
	end
end

%% do IFT and display
for ic = 1:ncoils
	ims(:,:,ic) = fftshift(ifftn(fftshift(d2d(:,:,ic))));
end
imsos = sqrt(sum(abs(ims).^2,3)); 
imagesc(imsos'); axis equal off; colormap gray;
title(sprintf('gradDelay = %.1f us, phsOffset = %.2f rad', gradDelay*1e6, phsOffset));

%% Estimate grad/acquisition delay and odd/even phase offset from image data (assumes that center strip is unaliased)
if 0
d2dOdd = d2d;
d2dOdd(:,2:2:end,:) = 0;
d2dEven = d2d;
d2dEven(:,1:2:end,:) = 0;
mask = imsos > 0.1*max(imsos(:));
for ic = 1:ncoils
	imsOdd(:,:,ic) = mask.*fftshift(ifftn(fftshift(d2dOdd(:,:,ic))));
	imsEven(:,:,ic) = mask.*fftshift(ifftn(fftshift(d2dEven(:,:,ic))));
end
pc2d = toppe.utils.phasecontrastMulticoil(imsEven,imsOdd);   % phase difference image (2D)
stripInds = (N/2-2):(N/2+1);                                 % unaliased center strip
pc1d = unwrap(mean(pc2d(:,(N/2-2):(N/2+2)), 2));
X = [ones(N,1) linspace(-pi,pi,N)'];
b = X\pc1d;
fprintf('Estimated shift: %.2f samples \nEstimated phase offset: %.2f radians \n', b(2), b(1));
end

return;

