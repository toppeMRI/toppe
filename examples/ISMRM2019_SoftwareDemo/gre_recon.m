function [ims imsos d]= gre_recon(pfile, readoutfile)
% Reconstruct 2D EPI data acquired with ISMRM2019 "live" demo
%
% Output:
%  ims           [nx ny ncoils]    
%  imsos         coil-combined (root-sum-of-squares) image
%  d             raw (k-space) data

addpath ~/gitlab/toppe/
%import toppe.*
%import toppe.utils.*

if ~exist('pfile','var')
	error('Usage: gre_recon(pfile, [readoutfile])');
end

if ~exist('readoutfile','var')
	readoutfile = 'readout.mod';
end

%% get readout file header
[~,gx,~,~,~,hdrints] = toppe.readmod(readoutfile);
ndat = size(gx,1);

%% load raw data
d = toppe.utils.loadpfile(pfile);   % int16, size [ndat ncoils nslices nechoes nviews] = [ndat ncoils 1 1 ny]
d = permute(d,[1 5 3 2 4]);         % [ndat ny nz ncoils nechoes]
d = double(d);
d = flipdim(d,1);        % data is stored in reverse order (for some reason)

%% get flat portion of readout
[rf,gx,gy,gz,desc,paramsint16,paramsfloat] = toppe.readmod(readoutfile);
nbeg = paramsint16(1);
nx = paramsint16(2);  % number of acquired data samples per TR, on plateau
decimation = round(125/paramsfloat(20));
d = d(nbeg:(nbeg+nx-1),:,:,:,:);     % [nx*125/oprbw ny nz ncoils nechoes]

%% Do IFT and display
for ic = 1:size(d,4)
	imstmp = fftshift(ifftn(fftshift(squeeze(d(:,:,1,ic)))));
   imstmp = imstmp(end/2+((-nx/decimation/2):(nx/decimation/2-1))+1,:,:);               % [nx ny nz]
	ims(:,:,ic) = imstmp;
end
imsos = sqrt(sum(abs(ims).^2,3)); 
im(imsos);

return;

