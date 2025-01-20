function [dat, header] = loadsafile(safile,varargin)
% function [dat, header] = loadsafile(safile,varargin)
%
% Load data from ScanArchive file.
%
% Input:
%   safile - ScanArchive filename string

% Output:
%   dat - raw k-space data from ScanArchive file
%   header - header info (currently just data size)
%
% Optional Inputs:
%   quiet - boolean, set to true to suppress command line output
%   acq_order - boolean, set to true to return the data in the order acq.
%   header_only - boolean, when true only returns header data (much faster)
%
% Default output dimensions of dat:
%  [ndat,ncoil,nslice,nechos,nview]
%
% Output dimensions of dat when 'acq_order' == true:
%  [ndat,ncoil,nframes]
% this will return the data in the order in which it was acquired



% (c) 2021 The Regents of the University of Michigan
% Melissa Haskell, mhask@umich.edu

import toppe.utils.*

%% Load input arguments
% Set defaults and parse varargin
arg.quiet        = false;
arg.acq_order    = false;
arg.header_only  = false;
arg.echo         = 'all'; % return all echos by default, but can be indexed
arg.version      = 'tv6';
arg = vararg_pair(arg, varargin);


%% make sure orchestra in on the path
% (note: not sure this is the best way to check this)
if exist('GERecon','file') ~= 3
    error('Error: loadsafile.m requires the GE orchestra toolbox.')
end


%% Load ScanArchive file
archive = GERecon('Archive.Load', safile);   


%% Load data size from header info

ndat = archive.DownloadData.rdb_hdr_rec.rdb_hdr_frame_size;
ncoil = archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).stop_rcv - ...
    archive.DownloadData.rdb_hdr_rec.rdb_hdr_dab(1).start_rcv + 1;
nsli = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nslices - 1; % "- 1"  for TOPPE slice issue
if strcmp(arg.version, 'tv7')
    nsli = nsli + 1;
end
necho = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nechoes;
nview = archive.DownloadData.rdb_hdr_rec.rdb_hdr_nframes;
nframe = archive.FrameCount;

if nframe ~= nsli*necho*nview
    warning('Warning: Check frame & data dimension sizes.'); 
end

% create header info (this could be expanded, put in dimension size for
% now)
header.dat_size = [ndat, ncoil, nsli, necho, nview];
if arg.header_only
    dat = []; 
    fprintf(1,'ndat = %d, nslices = %d, nechoes = %d, nviews = %d, ncoils = %d\n',...
        ndat, nsli, necho, nview, ncoil);
    return; 
end


%% Start loading in raw data

currentControl = GERecon('Archive.Next', archive);
kdata_tmp = currentControl.Data;  % size = [ndat ncoils]
ndat_cur_kdata = size(currentControl.Data, 1);
ncoil_cur_kdata = size(currentControl.Data, 2);

if ndat ~= ndat_cur_kdata, warning('Warning: Check data size.'); end
if ncoil ~= ncoil_cur_kdata, warning('Warning: Check data size.'); end

% init data matrix and write first frame to array
if arg.acq_order == true
    dat = zeros(ndat, ncoil, nframe);
    dat(:,:,1) = kdata_tmp;
else
    dat = zeros(ndat, ncoil, nsli, necho, nview);
    dat(:,:,currentControl.sliceNum,currentControl.echoNum+1,currentControl.viewNum) = kdata_tmp;
end


%% Load rest of frames

% loop through rest of frames
if ~arg.quiet
    fprintf(1,'ndat = %d, nslices = %d, nechoes = %d, nviews = %d, ncoils = %d\n',...
        ndat, nsli, necho, nview, ncoil);
    textprogressbar('Loading data: ');
end

if arg.acq_order == true
    for ii = 2:nframe
        if ~arg.quiet, textprogressbar(ii/(nview*nsli*necho)*100); end
        currentControl = GERecon('Archive.Next', archive);
        dat(:,:,ii) = currentControl.Data;   % size = [ndat ncoils nframes]
    end
else
    for ii = 2:nframe
        if ~arg.quiet, textprogressbar(ii/(nview*nsli*necho)*100); end
        currentControl = GERecon('Archive.Next', archive);
        kdata_tmp = currentControl.Data;
        iView = currentControl.viewNum;
        iSli = currentControl.sliceNum;
        iEcho = currentControl.echoNum+1;
        dat(:,:,iSli,iEcho,iView) = kdata_tmp;
    end
end
if ~arg.quiet, textprogressbar(' done.'); end

if ~strcmp(arg.echo,'all') 
    if arg.acq_order == true
        warning(['Cannot select echo indices when returning',...
            ' data in acquisition order, returning all echos.'])
    else
        dat = dat(:,:,:,arg.echo,:);
    end
end

end
