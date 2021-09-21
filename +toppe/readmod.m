function [rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = readmod(fname,showinfo)
% Read waveforms, ASCII description, and header arrays from TOPPE .mod file.
%
% function [rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = readmod(fname,showinfo)

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
%
% $Id: readmod.m,v 1.5 2018/11/13 18:07:29 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/readmod.m,v $

import toppe.*
import toppe.utils.*

if ~exist('showinfo','var')
	showinfo = false;
end

fid = fopen(fname, 'r', 'ieee-be');

% read ASCII description
asciisize = fread(fid, 1,         'int16');
desc      = fread(fid, asciisize, 'uchar'); 
desc = char(desc');

% read rest of header
hdr.ncoils = fread(fid, 1, 'int16');
hdr.res    = fread(fid, 1, 'int16');
hdr.npulses = fread(fid, 1, 'int16');
hdr.b1max  = fscanf(fid, 'b1max:  %f\n');
hdr.gmax   = fscanf(fid, 'gmax:   %f\n');

nparamsint16 = fread(fid, 1,            'int16');
paramsint16  = fread(fid, nparamsint16, 'int16');
nparamsfloat = fread(fid, 1,            'int16');
%paramsfloat  = fscanf(fid, '%f\n');                % this reads all float32 params at once
for n = 1:nparamsfloat
	paramsfloat(n)  = fscanf(fid, '%f\n', 1);
end

% get nChop (added in v4)
hdr.npre = paramsint16(1);  % number of discarded RF/ADC samples at start
hdr.rfres = paramsint16(2); % total number of RF/ADC samples

if showinfo
	fprintf(1, '\n%s', desc);
	fprintf(1, 'number of coils/channels:  %d\n', hdr.ncoils);
	fprintf(1, 'number points in waveform: %d\n', hdr.res);
	fprintf(1, 'number of waveforms:          %d\n', hdr.npulses);
	fprintf(1, 'data offset (bytes):       %d\n', ftell(fid));
	fprintf(1, '\n');
end

% read waveforms
rho = zeros(hdr.res,hdr.npulses,hdr.ncoils);
theta = rho;
gx = zeros(hdr.res,hdr.npulses);
gy = zeros(hdr.res,hdr.npulses);
gz = zeros(hdr.res,hdr.npulses);
for ip = 1:hdr.npulses
	for ic = 1:hdr.ncoils
		rho(:,ip,ic) = fread(fid, hdr.res, 'int16');
	end
	for ic = 1:hdr.ncoils
		theta(:,ip,ic) = fread(fid, hdr.res, 'int16');
	end
	gx(:,ip) = fread(fid, hdr.res, 'int16');
	gy(:,ip) = fread(fid, hdr.res, 'int16');
	gz(:,ip) = fread(fid, hdr.res, 'int16');
end

% convert back to physical units
max_pg_iamp = 2^15-2;                           % max instruction amplitude (max value of signed short)
rho   = rho*hdr.b1max/max_pg_iamp;     % Gauss
theta = theta*pi/max_pg_iamp;      % radians
gx    = gx*hdr.gmax/max_pg_iamp;       % Gauss/cm
gy    = gy*hdr.gmax/max_pg_iamp;
gz    = gz*hdr.gmax/max_pg_iamp;

% reshape output
rho   = reshape(rho,   hdr.res, hdr.npulses, hdr.ncoils);
theta = reshape(theta, hdr.res, hdr.npulses, hdr.ncoils);
gx    = reshape(gx,    hdr.res, hdr.npulses);
gy    = reshape(gy,    hdr.res, hdr.npulses);
gz    = reshape(gz,    hdr.res, hdr.npulses);
paramsint16 = paramsint16(4:end)';    % NB! Return only the user-defined ints passed to writemod.m
paramsfloat = paramsfloat';

rf = rho.*exp(1i*theta);

fclose(fid);



% EOF
