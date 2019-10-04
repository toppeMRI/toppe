function [dat, rdb_hdr] = loadpfile(pfile,echo,slicestart,sliceend,varargin)
% function [dat, rdb_hdr] = loadpfile(pfile,[echo,slicestart])
%
% Load data for one echo (or all) from Pfile, EXCEPT dabslice=0 slot (which can contain corrupt data).
%
% Input options:
%  echo          only get data for this echo (default: load all echoes)
%  slicestart    get data starting from this slice index (1:N slices, default: 2)
%  sliceend      "" except ending slice (default: N)
%
% Output dimensions of dat:
%  [ndat,ncoil,nslice,nechos,nview]
%
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

%% Load input arguments
% Set defaults and parse varargin
arg.quiet        = false;
arg = vararg_pair(arg, varargin);

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

%% Determine which echoes to load
if exist('echo','var') && ~isempty(echo)
	ECHOES = echo;
else
	ECHOES = 1:nechoes;
end
if nargin < 3
	slicestart = 2;
end

if nargin < 4
    sliceend = nslices;
end

%% Calculate size of data chunks. See pfilestruct.jpg, and rhrawsize calculation in .e file.
echores  = ndat*(nviews+1);                 % number of data points per 'echo' loaddab slot. Includes baseline (0) view.
sliceres = nechoes*echores;                 % number of data points per 'slice'
coilres  = nslices*sliceres;                % number of data points per receive coil

pfilesize = rdb_hdr.off_data + 2*ptsize*ncoils*nslices*nechoes*(nviews+1)*ndat;   % this should match the Pfile size exactly
pfilename=dir(pfile);

% File size check for integrity
if pfilesize ~= pfilename.bytes
    warning('Expected %0.1fMB file but read in %0.1fMB file.\n',pfilesize/1e6,pfilename.bytes/1e6)
    fprintf('Press enter to continue anyway...');
    input('');
end

% Check if second view is empty
% This happens when there's only one view, but the scanner sets nviews = 2
if nviews == 2
    fseek(fid, rdb_hdr.off_data+(2*ptsize*(sliceres + 2*ndat)), 'bof'); % Seek to view slice 2, echo 1, view 2, coil 1
    view2tmp = fread(fid, 2*ndat, 'int16=>int16');
    if all(view2tmp==0)
        nviews = 1; % Set nviews to 1 so we don't read in all the empty data
        if ~arg.quiet
            warning('off','backtrace'); % Turn off backtrace lines for cleaner output
            warning('View 2 appears to be empty, only loading view 1.');
        end
    end
end

% Try to check how much free memory we have (linux implementation) and warn if
% we are going to exceed it. If we exceed free RAM we may cause a system
% freeze due to attempted mem swap

try
    % Check predicted size of output vs free memory
    memneeded = 16*ndat*ncoils*(sliceend-slicestart+1)*numel(ECHOES)*nviews; %16 bytes per complex value;
    
    % Pull available memory using linux system command
    [~,out]=system('cat /proc/meminfo | grep "MemAvailable:"');
    memfree=sscanf(out,'MemAvailable: %f'); % Available memory in kB
    mempercentuse = (memneeded/1000) / memfree;
    if mempercentuse > 1
        warning('Loading data (%0.1fGB) may exceed available RAM. This could end up freezing your computer!!',memneeded/1e9);
        fprintf('Press enter to continue anyway...');
        input('');
    elseif mempercentuse > 0.9 % Check if we are going to use 90% of memory and warn
        warning('Loading data (%0.1fGB) is going to use %0.1f%% of your available RAM. Proceed with caution!',memneeded/1e9,100*mempercentuse);
    end
catch
end

if ~arg.quiet
fprintf(1,'ndat = %d, nslices = %d, nechoes = %d, nviews = %d, ncoils = %d\n', ndat, nslices, nechoes, nviews, ncoils);
end

%% Read data from file
datr = int16(zeros(ndat,ncoils,sliceend-slicestart+1,numel(ECHOES),nviews));
dati = datr;
if ~arg.quiet; textprogressbar('Loading data: '); end
for icoil = 1:ncoils
    if ~arg.quiet; textprogressbar(icoil/ncoils*100); end
    for islice = slicestart:sliceend   % skip first slice (sometimes contains corrupted data)
        for iecho = 1:length(ECHOES) % Load every element in ECHOES
				echo = ECHOES(iecho);
            for iview = 1:nviews
                offsetres = (icoil-1)*coilres + (islice-1)*sliceres + (echo-1)*echores + iview*ndat;
                offsetbytes = 2*ptsize*offsetres;
                fseek(fid, rdb_hdr.off_data+offsetbytes, 'bof');
                dtmp = fread(fid, 2*ndat, 'int16=>int16');
                sliceind = islice-slicestart+1;
                datr(:,icoil,sliceind,iecho,iview) = dtmp(1:2:end); %Real data
                dati(:,icoil,sliceind,iecho,iview) = dtmp(2:2:end); %Imag
            end
        end
    end
end

%Combine real+imag in one step. This ends up using 2 times as much memory
%since the complex data exists twice (datr+dati and dat) but is faster
%since complex function is vectorized.

dat = complex(datr,dati); % Combine data in one step
clearvars datr dati % Free up some memory
dat = double(dat);  % Convert to double in place
fclose(fid);
if ~arg.quiet; textprogressbar(' done.'); end
return
