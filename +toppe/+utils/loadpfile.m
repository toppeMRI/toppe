function [dat, rdb_hdr] = loadpfile(pfile,echo)
% function [dat, rdb_hdr] = loadpfile(pfile,[echo])
%
% Load data for one echo (or all) from Pfile, EXCEPT dabslice=0 slot (which can contain corrupt data).

% This file is part of the TOPPE development environment for platform-independent MR pulse programming.
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

import toppe.utils.*

%% Loadpfile code
fid = fopen(pfile,'r','l');
ver = fread(fid,1,'float32');
str = num2str(ver);
rdbm_rev = str2double(str);
fseek(fid,0,'bof');                 % NB!
rdb_hdr = read_rdb_hdr(fid,rdbm_rev);


%% Header parameters
ndat    = rdb_hdr.frame_size;
nslices = rdb_hdr.nslices;
ptsize  = rdb_hdr.point_size;                    % Either 2 (data stored in short int format, int16) or 4 (extended precision)
nechoes = rdb_hdr.nechoes;
nviews  = rdb_hdr.nframes;
ncoils  = rdb_hdr.dab(2)-rdb_hdr.dab(1)+1;

%% Determine which echoes to load in
if exist('echo','var')
	ECHOES = echo;
else
	ECHOES = 1:nechoes;
end

%% Calculate size of data chunks. See pfilestruct.jpg, and rhrawsize calculation in .e file.
echores  = ndat*(nviews+1);                 % number of data points per 'echo' loaddab slot. Includes baseline (0) view.
sliceres = nechoes*echores;                 % number of data points per 'slice'
coilres  = nslices*sliceres;                % number of data points per receive coil

pfilesize = rdb_hdr.off_data + 2*ptsize*ncoils*nslices*nechoes*(nviews+1)*ndat;   % this should match the Pfile size exactly
pfilename=dir(pfile);

if pfilesize ~= pfilename.bytes
    warning('Expected %0.1fMB file but read in %0.1fMB file.\n',pfilesize,pfilename.bytes/1e6)
    fprintf('Press enter to continue anyway...');
    input('');
end

fprintf(1,'\nndat = %d, nslices = %d, nechoes = %d, nviews = %d, ncoils = %d\n', ndat, nslices, nechoes, nviews, ncoils);

%% Read data from file
datr = int16(zeros(ndat,ncoils,nslices-1,numel(ECHOES),nviews));
dati = datr;
textprogressbar('Loading data: ');
for icoil = 1:ncoils
    textprogressbar(icoil/ncoils*100);
    for islice = 2:nslices   % skip first slice (sometimes contains corrupted data)
        for iecho = ECHOES % Load every element in ECHOES
            for iview = 1:nviews
                offsetres = (icoil-1)*coilres + (islice-1)*sliceres + (iecho-1)*echores + iview*ndat;
                offsetbytes = 2*ptsize*offsetres;
                fseek(fid, rdb_hdr.off_data+offsetbytes, 'bof');
                dtmp = fread(fid, 2*ndat, 'int16=>int16');
                datr(:,icoil,islice-1,iecho,iview) = dtmp(1:2:end); %Real data
                dati(:,icoil,islice-1,iecho,iview) = dtmp(2:2:end); %Imag
            end
        end
    end
end
%Combine real+imag in one step. This ends up using 2 times as much memory
%since the complex data exists twice (datr+dati and dat) but is faster
%since complex function is vectorized. May cause issues if your computer is
%RAM limited...

dat = complex(datr,dati); % Combine data in one step
clearvars datr dati % Free up some memory
dat = double(dat);  % Convert to double in place
fclose(fid);
textprogressbar(' done.');
return
