function [dat, pfilesize] = loaddat_ge(fid,rdb_hdr,SLICES,ECHOES,VIEWS,doflip,COILS)
% function [dat, pfilesize] = loaddat_ge(fid,rdb_hdr,slice,echo,view, [doflip,COILS])
%
% Load k-space data from a GE raw data file (Pfile).
% 
% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
%
% OUTPUT:
%  dat:      [ndat ncoils] array containing data for one view (all coils).
%
% INPUTS:
%  fid:      Pfile id, e.g., created with 
%               fid = fopen(pfile,'r','l');
%  rdb_hdr:  Pfile header, created with
%               fseek(fid,0,'bof');
%               rdb_hdr = read_rdb_hdr(fid,rdbm_rev);   % the fMRI lab at UM is currently at rdbm_rev = 24
%            To get the Pfile version, do:
%               fseek(fid,0,'bof');
%               ver = fread(fid,1,'float32');
%               str = num2str(ver);
%               rdbm_rev = str2double(str);
%  slice:    'slice' index (in loaddab), starting at 0
%  echo:     'echo' index (in loaddab), starting at 0
%  view:     'view' index (in loaddab), starting at 0
%  doflip:   1: reverse data order. 0: don't. (optional. default = 0)
%  COILS:    Return data for these coils (optional. default = all coils)
%
% TOPPE is free software: you can redistribute it and/or modify
% it under the terms of the GNU Library General Public License as published by
% the Free Software Foundation version 2.0 of the License.
%
% TOPPE is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Library General Public License for more details.
%
% You should have received a copy of the GNU Library General Public License
% along with TOPPE. If not, see <http://www.gnu.org/licenses/old-licenses/lgpl-2.0.html>.
% 
% (c) 2016 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu
%
% $Id: loaddat_ge.m,v 1.2 2018/10/24 19:33:34 jfnielse Exp $

ncoils  = rdb_hdr.dab(2)-rdb_hdr.dab(1)+1;
ptsize  = rdb_hdr.point_size;                    % Either 2 (data stored in short int format, int16) or 4 (extended precision)
nviews  = rdb_hdr.nframes;
nechoes = rdb_hdr.nechoes;
nslices = rdb_hdr.nslices;

% Calculate size of data chunks. See pfilestruct.jpg, and rhrawsize calculation in .e file.
ndat     = rdb_hdr.frame_size;              % number of data points per view
echores  = ndat*(nviews+1);                 % number of data points per 'echo' loaddab slot. Includes baseline (0) view.
sliceres = nechoes*echores;                 % number of data points per 'slice'
coilres  = nslices*sliceres;                % number of data points per receive coil

pfilesize = rdb_hdr.off_data + 2*ptsize*ncoils*nslices*nechoes*(nviews+1)*ndat;   % this should match the Pfile size exactly

% read data
dat = zeros([ndat ncoils]);
if ~exist('COILS','var') 
	COILS = 1:ncoils;
end
dat = int16(zeros(ndat,numel(COILS),numel(SLICES),numel(ECHOES),numel(VIEWS)));
for ic = 1:numel(COILS)
	coil = COILS(ic);
	for isl = 1:numel(SLICES)
		slice = SLICES(isl);
		for ie = 1:numel(ECHOES)
			echo = ECHOES(ie);
			for iv = 1:numel(VIEWS)
				view = VIEWS(iv);
				offsetres = (coil-1)*coilres + slice*sliceres + echo*echores + view*ndat;
				offsetbytes = 2*ptsize*offsetres;
				fseek(fid, rdb_hdr.off_data+offsetbytes, 'bof');
				if ptsize == 2
					dtmp = fread(fid, 2*ndat, 'int16');
				else
					dtmp = fread(fid, 2*ndat, 'int32');   % NB! Has not been tested
				end
				dtmp = complex(dtmp(1:2:end), dtmp(2:2:end));
				dat(:,ic,isl,ie,iv) = int16(dtmp);
			end
		end
	end
end
%ftell(fid)

dat = squeeze(dat);

if exist('doflip','var') 
	if doflip == 1
		dat = flipdim(dat,1);
	end
end

return;
