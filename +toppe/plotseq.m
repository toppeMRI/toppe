function [rf, gx, gy, gz, rf1, gx1, gy1, gz1, tdelay] = plotseq(nstart, nstop, varargin)
% Display pulse sequence, as specified in modules.txt, scanloop.txt, and timing.txt
%
% function [rf, gx, gy, gz, rf1, gx1, gy1, gz1, tdelay] = plotseq(nstart, nstop, varargin)
%
% Inputs:
%   nstart,nstop       first and last startseq calls (as specified in scanloop.txt)
% Options:
%   loopFile           default: 'scanloop.txt'
%   loopArr            scan loop array (see readloop.m). Default: read from loopFile.
%   moduleListFile     default: 'modules.txt'
%   mods               Structure containing .mod file contents (see ./utils/readModulesListFile.m). 
%                      Use in playseq.m to speed up display.
%                      Default: get values from moduleListFile.
%   doDisplay          true (default) or false
%   system             struct specifying hardware system info, see systemspecs.m
%
% Outputs:
%   rf               Complex RF waveform (Gauss)
%   gx,gy,gz         Gauss/cm
%   rf1/gx1/...      Values from most recent startseq() call (used in, e.g., ge2seq.m and playseq.m)
%   tdelay           Delay after end of module waveform. 
%                    Determined by duration in modules.txt AND by textra (column 14) in scanloop.txt   
%                    Used in ge2seq.m to set delay block, and in playseq.m.
%
% $Id: plotseq.m,v 1.6 2018/11/05 00:12:32 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/plotseq.m,v $

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
% (c) 2016-2018 The Regents of the University of Michigan
% Jon-Fredrik Nielsen, jfnielse@umich.edu

%% parse inputs

% Default values 
arg.loopArr         = [];
arg.loopFile        = 'scanloop.txt';
arg.mods            = [];
arg.moduleListFile = 'modules.txt';
arg.doDisplay       = true;
arg.system          = toppe.systemspecs();  % Accept default timing (includes EPIC-related time gaps)

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% read scan files as needed

% scanloop array
if isempty(arg.loopArr)
	loopArr = toppe.utils.tryread(@toppe.readloop, arg.loopFile);
else
	loopArr = arg.loopArr;
end

% module waveforms
if isempty(arg.mods)
	cores = toppe.utils.tryread(@toppe.readmodulelistfile, arg.moduleListFile);
else
	cores = arg.mods;
end

%% timing CVs
c = struct2cell(arg.system.toppe);
TPARAMS = cell2mat(c(2:end));
[start_core myrfdel daqdel timetrwait timessi] = deal(TPARAMS(1), TPARAMS(2), TPARAMS(3), TPARAMS(4), TPARAMS(5));

%% build sequence. each sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
dt = 4;  % us
max_pg_iamp = 2^15-2;
for it = nstart:nstop
	ic = loopArr(it,1);   % core id 
	if cores{ic}.hasRF
		ia_rf = loopArr(it,2);
	else
		ia_rf = 0;
	end
	ia_th = loopArr(it,3);
	ia_gx = loopArr(it,4);
	ia_gy = loopArr(it,5);
	ia_gz = loopArr(it,6);

	if cores{ic}.hasRF
		coredel = myrfdel;
	elseif cores{ic}.hasDAQ
		coredel = daqdel;
	else
		coredel = 0;
	end

	tmin = start_core + coredel + cores{ic}.wavdur + timetrwait + timessi;   % mimimum core duration (us). 
	tdelay = max(cores{ic}.dur - tmin, 0);                                   % silence at end of core
	tminwait = 12;   % (us) min length of wait pulse.
	if size(loopArr,2)>13
		tdelay = tdelay + max(loopArr(it,14),tminwait);    % waitcore duration (see toppev2.e)
	end

	waveform = loopArr(it,16);

	% get gradients and apply in-plane (xy) rotation
	gxit = cores{ic}.gx(:,waveform);
	gyit = cores{ic}.gy(:,waveform);
	gzit = cores{ic}.gz(:,waveform);
	iphi = loopArr(it,11);
	phi = iphi/max_pg_iamp*pi;    % rad, [-pi pi]
	Gxy = [cos(phi) -sin(phi); sin(phi) cos(phi)]*[gxit(:)'; gyit(:)'];
	gxit = Gxy(1,:)';
	gyit = Gxy(2,:)';
	
	rho1 = [zeros(round((start_core+coredel)/dt),1); ia_rf/max_pg_iamp*  abs(cores{ic}.rf(:,waveform));  zeros(round((timetrwait+timessi)/dt),1)];
	th1  = [zeros(round((start_core+coredel)/dt),1); ia_th/max_pg_iamp*angle(cores{ic}.rf(:,waveform));  zeros(round((timetrwait+timessi)/dt),1)];
	gx1  = [zeros(round((start_core)/dt),1);         ia_gx/max_pg_iamp*gxit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gy1  = [zeros(round((start_core)/dt),1);         ia_gy/max_pg_iamp*gyit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];
	gz1  = [zeros(round((start_core)/dt),1);         ia_gz/max_pg_iamp*gzit(:); zeros(round((timetrwait+timessi+coredel)/dt),1)];

	% apply RF phase offset
	if cores{ic}.hasRF
		th1 = th1 + loopArr(it,12)/max_pg_iamp*pi;
		th1 = angle(exp(1i*th1));   % wrap to [-pi pi] range
	end

	rho = [rho; rho1; zeros(round(tdelay/dt),1)];
	th  = [th;  th1;  zeros(round(tdelay/dt),1)];
	gx  = [gx;  gx1;  zeros(round(tdelay/dt),1)];
	gy  = [gy;  gy1;  zeros(round(tdelay/dt),1)];
	gz  = [gz;  gz1;  zeros(round(tdelay/dt),1)];

	%fprintf(1, 'it %d: tmin = %.3f ms, rf t = %.3f ms, grad t = %.3f ms\n', it, tmin/1000, numel(rho)*dt*1e-3, numel(gx)*dt*1e-3);
end

rf = rho.*exp(1i*th);
rf1 = rho1.*exp(1i*th1);     % waveforms in last module, without the delay after it (if any)

% plot
if arg.doDisplay
	T = (0:(numel(rho)-1))*dt/1000; % msec
	gmax = 5;  % Gauss/cm
    srho = max(1.1*max(abs(rho(:))),0.05);
	subplot(511); plot(T, rho); ylabel('rho');   axis([T(1) 1.01*T(end) -srho srho]);
	subplot(512); plot(T, th);  ylabel('theta'); axis([T(1) 1.01*T(end) -1.3*pi 1.3*pi]);
	subplot(513); plot(T, gx);  ylabel('gx'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	%gmax = 1;  % Gauss/cm
	subplot(514); plot(T, gy);  ylabel('gy'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	subplot(515); plot(T, gz);  ylabel('gz'); axis([T(1) 1.01*T(end) -1.05*gmax 1.05*gmax]);;
	xlabel('msec');
end

return;

% EOF
