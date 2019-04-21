function [ims imsos d]= epi_recon(pfile, readoutfile)
% Reconstruct 2D EPI data acquired with ISMRM2019 "live" demo
%
% Output:
%  ims           [nx ny ncoils]    
%  imsos         coil-combined (root-sum-of-squares) image
%  d             raw (k-space) data

test = true;

import toppe.*
import toppe.utils.*

if ~exist('readoutfile','var')
	readoutfile = 'readout.mod';
end

% get readout file header
[~,gx,~,~,~,hdrints] = toppe.readmod(readoutfile);
ndat = size(gx,1);
N    = hdrints(1);       % image size
nes  = hdrints(3);       % echo spacing (number of 4us samples)
npre = hdrints(4);       % number of samples before start of readout plateau of first echo
nshots = size(gx,2);

% load raw data
if ~test
d = loadpfile(pfile);               % int16, size [ndat ncoils nslices nechoes nviews] = [ndat ncoils 1 1 nshots]
d = permute(d,[1 5 2 3 4]);         % [ndat nshots ncoils].
d = double(d);
d = flipdim(d,1);        % data is stored in reverse order for some reason
[ndat nshots ncoils] = size(d);
else
	ncoils = 1;
end
	

% sort data into 2D NxN matrix
d2d = zeros(N,N,ncoils);
etl = N/nshots;    % echo-train length
iy = 1;
for ic = 1:ncoils
	for ii = 1:1 %nshots
		for jj = 1:etl
			istart = npre + (jj-1)*nes + 1;
			d2d(:,iy,ic) = d(istart:(istart+N-1), ii, ic);
		end
	end
end

return

for coil = 1:size(d,4)
	im(:,:,coil) = fftshift(ifftn(fftshift(D)));
	imstmp = imstmp(end/2+((-nx/decimation/2):(nx/decimation/2-1))+1,:,:);               % [nx ny nz]
	ims(:,:,:,coil) = imstmp;
end

fprintf('\n');

%ims = flipdim(ims,1);

% display root sum-of-squares image
imsos = sqrt(sum(abs(ims).^2,4)); 
%figure; im('blue0',imsos,[0 1.3]);
if exist('clim','var')
	if dodisplay
		%im(permute(imsos,[2 1 3]),clim);
		im(imsos,clim);
	end
else
	if dodisplay
		%im(permute(imsos,[2,1,3]));
		im(imsos);
	end
end

return;

