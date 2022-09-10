function [ims imsos d]= recon3dft(pfile,varargin) 
% Recon 3D spin-warp image acquired with TOPPE.
%
% function [ims imsos d]= recon3dft(pfile,varargin) 
%
% Beginning of duration of readout plateau (in samples) is stored in readout.mod header, 
% See, e.g., makeepi.m
%
% Input options:
%  type             '2d' or '3d'
%  echo             echo to recon. Default: 1
%  readoutFile      default: 'readout.mod'
%  dokzft           default: true
%  zpad             default: [1 1]
%  dodisplay        default: false
%  flipFid          default: true. Flip image along readout (fid) dimension to match scanner display
%  alignWithUCS     default: false. Flip images along 2nd and 3rd dimension so it appears correctly
%                   in 'im' MIRT viewer. UCS = 'universal coordinate system'
%
% Output:
%  ims           [nx ny nz ncoils]    
%  imsos         [nx ny nz]   root-sum-of-squares coil-combined image

import toppe.*
import toppe.utils.*

%echo,readoutFile,dokzft,zpad,dodisplay,clim)
arg.type = '3d';
arg.echo = 1;
arg.readoutFile = 'readout.mod';
arg.dokzft = 'true';
arg.zpad = [1 1];
arg.dodisplay = false;
arg.flipFid = true;
arg.alignWithUCS = false;

arg = toppe.utils.vararg_pair(arg, varargin);

if arg.alignWithUCS
    if ~arg.flipFid
        warning('flipFid set to true to align with UCS');
    end
    arg.flipFid = true; 
end

echo = arg.echo;
readoutFile = arg.readoutFile;
zpad = arg.zpad;
dodisplay = arg.dodisplay;

% load raw data
d = loadpfile(pfile,echo);   % int16. [ndat ncoils nslices nechoes nviews] = [ndat ncoils nz 2 ny]
d = permute(d,[1 5 3 2 4]);         % [ndat ny nz ncoils nechoes]
d = double(d);

d = flipdim(d,1);        % data is stored in reverse order for some reason

%if(mod(size(d,3),2))
%	d = d(:,:,2:end,:,:);  % throw away dabslice = 0
%end

% get flat portion of readout
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = readmod(readoutFile);
nramp = 0; %15;  % see writemod.m
nbeg = paramsint16(1) - hdr.npre + nramp;
nx = paramsint16(2);  % number of acquired data samples per TR
decimation = round(125/paramsfloat(20));
d = d(nbeg:(nbeg+nx-1),:,:,:,:);     % [nx*125/oprbw ny nz ncoils nechoes]

% zero-pad in z
if zpad(2) > 1
	[ndat ny nz ncoils nechoes] = size(d);
	d2 = zeros([ndat ny round(nz*zpad(2)) ncoils]);
	d2(:,:,(end/2-nz/2):(end/2+nz/2-1),:,:) = d;
	d = d2; clear d2;
end

% recon 
for coil = 1:size(d,4)
	fprintf(1,'\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\brecon coil %d', coil);
	imstmp = ift3(d(:,:,:,coil), 'type', arg.type);
	imstmp = imstmp(end/2+((-nx/decimation/2):(nx/decimation/2-1))+1,:,:);               % [nx ny nz]
	if zpad(1) > 1   % zero-pad (interpolate) in xy
		dtmp = fft3(imstmp);
		[nxtmp nytmp nztmp] = size(dtmp);
		dtmp2 = zeros([round(nxtmp*zpad(1)) round(nytmp*zpad(1)) nztmp]);
		dtmp2((end/2-nxtmp/2):(end/2+nxtmp/2-1),(end/2-nytmp/2):(end/2+nytmp/2-1),:) = dtmp;
		dtmp = dtmp2; clear dtmp2;
		imstmp = ift3(dtmp);
	end
	ims(:,:,:,coil) = imstmp;
end

fprintf('\n');

% Flip along L/R to match the scanner host display.
% In Axial view on console: 'R' is on left; 'A' is top
% In Sagittal view: 'A' is on left; 'S' is on top
% In Coronal view: 'R' is on left; 'S' is on top
if arg.flipFid
    ims = flipdim(ims,1);
    warning('First image dimension flipped to match host display.');
end

if arg.alignWithUCS
    ims = flipdim(ims,2);
    ims = flipdim(ims,3);
end

% root sum-of-squares image
imsos = sqrt(sum(abs(ims).^2,4)); 

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


function im = sub_ift3(D,do3dfft)
%
%	function im = sub_ift3(dat)
%
%	Centered inverse 3DFT of a 3D data matrix.
% 
% $Id: recon3dft.m,v 1.8 2018/11/13 18:41:36 jfnielse Exp $

if ~exist('do3dfft','var')
	do3dfft = true;
end

if do3dfft
	im = fftshift(ifftn(fftshift(D)));
else
	% don't do fft in 3rd dimension
	for k = 1:size(D,3)
		im(:,:,k) = fftshift(ifftn(fftshift(D(:,:,k))));
	end
end
    

return;
