function [ims imsos d]= epi_recon(pfile, readoutfile)
% Reconstruct 2D EPI data acquired with ISMRM2019 "live" demo
%
% Output:
%  ims           [nx ny ncoils]    
%  imsos         coil-combined (root-sum-of-squares) image
%  d             raw (k-space) data

addpath ~/gitlab/toppe/
%import toppe.*
%import toppe.utils.*

if ~exist('readoutfile','var')
	readoutfile = 'readout.mod';
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

%% apply gradient/acquisition delay
d = circshift(d, 1);
%dup = interpft(d, 5*ndat
%for ii = 1:nshots
%	for ic = 1:ncoils
%		dtmp = 

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
im(imsos);
%keyboard

return;

