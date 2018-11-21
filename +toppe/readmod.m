function [rf,gx,gy,gz,desc,paramsint16,paramsfloat] = readmod(fname,showinfo)
% Read waveforms, ASCII description, and header arrays from TOPPE .mod file.
%
% function [rf,gx,gy,gz,desc,paramsint16,paramsfloat] = readmod(fname,showinfo)

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
ncoils = fread(fid, 1, 'int16');
res    = fread(fid, 1, 'int16');
npulses = fread(fid, 1, 'int16');
b1max  = fscanf(fid, 'b1max:  %f\n');
gmax   = fscanf(fid, 'gmax:   %f\n');

nparamsint16 = fread(fid, 1,            'int16');
paramsint16  = fread(fid, nparamsint16, 'int16');
nparamsfloat = fread(fid, 1,            'int16');
%paramsfloat  = fscanf(fid, '%f\n');                % this reads all float32 params at once
for n = 1:nparamsfloat
	paramsfloat(n)  = fscanf(fid, '%f\n', 1);
end

if showinfo
	fprintf(1, '\n%s', desc);
	fprintf(1, 'number of coils/channels:  %d\n', ncoils);
	fprintf(1, 'number points in waveform: %d\n', res);
	fprintf(1, 'number of waveforms:          %d\n', npulses);
	fprintf(1, 'data offset (bytes):       %d\n', ftell(fid));
	fprintf(1, '\n');
end

% read waveforms
rho = zeros(res,npulses,ncoils);
theta = rho;
gx = zeros(res,npulses);
gy = zeros(res,npulses);
gz = zeros(res,npulses);
for ip = 1:npulses
	for ic = 1:ncoils
		rho(:,ip,ic) = fread(fid, res, 'int16');
	end
	for ic = 1:ncoils
		theta(:,ip,ic) = fread(fid, res, 'int16');
	end
	gx(:,ip) = fread(fid, res, 'int16');
	gy(:,ip) = fread(fid, res, 'int16');
	gz(:,ip) = fread(fid, res, 'int16');
end

% convert back to physical units
max_pg_iamp = 2^15-2;                           % max instruction amplitude (max value of signed short)
rho   = rho*b1max/max_pg_iamp;     % Gauss
theta = theta*pi/max_pg_iamp;      % radians
gx    = gx*gmax/max_pg_iamp;       % Gauss/cm
gy    = gy*gmax/max_pg_iamp;
gz    = gz*gmax/max_pg_iamp;

% reshape output
rho   = reshape(rho,   res, npulses, ncoils);
theta = reshape(theta, res, npulses, ncoils);
gx    = reshape(gx,    res, npulses);
gy    = reshape(gy,    res, npulses);
gz    = reshape(gz,    res, npulses);
paramsint16 = paramsint16(3:end)';             % NB! Return only the user-defined ints passed to writemod.m
paramsfloat = paramsfloat';

rf = rho.*exp(1i*theta);

fclose(fid);



% EOF
