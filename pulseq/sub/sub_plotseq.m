function [rf,gx,gy,gz] = sub_plotseq(moduleArr, loopStructArr, nstart, nstop, varargin)
% Plot sequence contained in moduleArr and loopStructArr.
% Skips delay blocks (for now)
% See also sub_playseq.m

%% parse inputs
arg.nTRskip = 0;
arg.system    = toppe.systemspecs();  % Accept default timing (includes EPIC-related time gaps)
arg.gradMode  = 'amplitude';          % 'amplitude' or 'slew'
arg.doDisplay = true;
arg.units     = 'Gauss';              % or 'Hz', like Pulseq uses

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

%% timing CVs
c = struct2cell(arg.system.toppe);
TPARAMS = cell2mat(c(2:8));
[start_core_rf start_core_daq start_core_grad myrfdel daqdel timetrwait timessi] = ...
	deal(TPARAMS(1), TPARAMS(2), TPARAMS(3), TPARAMS(4), TPARAMS(5), TPARAMS(6), TPARAMS(7));

%% build sequence. each sample is 4us
rf = []; gx = []; gy = []; gz = [];
dt = 4;  % us
for it = nstart:nstop
	loopStruct = loopStructArr(it);

	if isempty(loopStruct.mod)
		% skip delay blocks. TODO: add
		continue;
	end

	module = moduleArr(loopStruct.mod);

	if module.hasRF
		coredel = myrfdel;
		start_core = start_core_rf;
	elseif module.hasADC
		coredel = daqdel;
		start_core = start_core_daq;
	else
		coredel = 0;
		start_core = start_core_grad;
	end

	tmin = start_core + coredel + dt*module.nt + timetrwait + timessi;   % mimimum core duration (us).
	tminwait = 12;   % (us) min length of wait pulse.
	tdelay = max(round(loopStruct.textra*1e6),tminwait);    % us

	wavnum = loopStruct.wavnum;

	if ~isempty(module.rf)
		% construct rf waveform with phase offset
		rf1 = loopStruct.rfamp*module.rf(:,wavnum);
		rf1 = rf1.*exp(1i*loopStruct.rfphs);
	else
		rf1 = zeros(module.nt,1);
	end

	for ax = {'gx','gy','gz'}
		ax = ax{1};
		eval(sprintf('ie = isempty(module.%s);', ax));
		if ~ie
			eval(sprintf('%s1 = loopStruct.%samp*module.%s(:,wavnum);', ax, ax, ax));
		else
			eval(sprintf('%s1 = zeros(module.nt,1);', ax));
		end
	end

	rf = [rf; zeros(round((start_core+coredel)/dt),1); rf1; zeros(round((tdelay+timetrwait+timessi)/dt),1)];
	gx = [gx; zeros(round((start_core)/dt),1); gx1; zeros(round((coredel+tdelay+timetrwait+timessi)/dt),1)];
	gy = [gy; zeros(round((start_core)/dt),1); gy1; zeros(round((coredel+tdelay+timetrwait+timessi)/dt),1)];
	gz = [gz; zeros(round((start_core)/dt),1); gz1; zeros(round((coredel+tdelay+timetrwait+timessi)/dt),1)];
end

%% Convert units
if strcmp(arg.units, 'Hz')
	rfscale = arg.system.gamma;                 % convert from Gauss to Hz
	gscale = 100*arg.system.gamma*1e-3;         % convert from Gauss/cm to kHz/m
	rfunit = 'Hz';
	gunit = 'kHz/m';
else
	rfscale = 1;
	gscale = 1;
	rfunit = 'Gauss';
	gunit = 'Gauss/cm';
end
rf = rfscale*rf;
gx = gscale*gx;
gy = gscale*gy;
gz = gscale*gz;

%% plot
if arg.doDisplay
   T = (0:(numel(rf)-1))*dt/1000; % msec
   Tend = 1.01*T(end);

   rho = abs(rf);
   th = angle(rf);

   % rf
   srho = max(1.1*max(abs(rho(:))),0.05);
   lw = 1.5;
   subplot(511); plot(T, rho, 'LineWidth', lw); ylabel(sprintf('rho (%s)',rfunit)); axis([T(1) Tend -srho srho]);
   title(sprintf('[%d,%d]', nstart, nstop));
   subplot(512); plot(T, th, 'LineWidth', lw);  ylabel('theta (rad)'); axis([T(1) Tend -1.3*pi 1.3*pi]);

	% gradients/slew
	if strcmp(arg.gradMode, 'slew')
		T = T(2:end);
		gx = diff(gx)/dt*1e3;   % G/cm/ms
		gy = diff(gy)/dt*1e3;
		gz = diff(gz)/dt*1e3;
	end
	gmax = max(abs([gx(:); gy(:); gz(:)]));
	subplot(513); plot(T, gx, 'LineWidth', lw);  ylabel(sprintf('gx (%s)',gunit)); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
	subplot(514); plot(T, gy, 'LineWidth', lw);  ylabel(sprintf('gy (%s)',gunit)); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
	subplot(515); plot(T, gz, 'LineWidth', lw);  ylabel(sprintf('gz (%s)',gunit)); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
	xlabel('msec');
end

return;

