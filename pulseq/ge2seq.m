function seq = ge2seq(toppeTarFile, varargin)
% function seq = ge2seq(toppeTarFile, varargin)
%
% TOPPE to Pulseq file conversion.
% Inputs:
%  toppeTarFile       TOPPE .tgz archive file containing the following:
%                     *.mod:            One .modfile corresponds to one "unique" block (see below)
%                     modules.txt       List of .mod files, and flags indicating whether the .wav file is RF/ADC/(gradients only)
%                     scanloop.txt      Sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)
% Options:
%  seqfile            Output .seq file name
%  moduleListFile     Text file listing all .mod files. Default: 'modules.txt'.
%                     The .mod files listed must exist in the Matlab path.
%  loopFile           Text file specifying the MR scan loop. Default: 'scanloop.txt'
%  system             struct specifying Siemens scanner hardware limits, see mr.opts
%  systemGE           struct specifying GE system specs, see +toppe/systemspecs.m
%
% Examples:
%  >> ge2seq('TOPPEseq.tar');
%
%  >> lims = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
%                    'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
%                    'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  
%  >> sys = systemspecs('maxSlew', 130, 'slewUnit', 'T/m/s');
%  >> ge2seq('TOPPEseq.tar', 'system', lims, 'systemGE', sys);
%


%% Parse inputs and set system values
% defaults
arg.seqfile        = 'out.seq';
arg.debug          = false;
arg.debugAdc       = false;
arg.moduleListFile = 'modules.txt';
arg.loopFile       = 'scanloop.txt';

arg.system = mr.opts('MaxGrad', 32, 'GradUnit', 'mT/m',...
                     'MaxSlew', 130, 'SlewUnit', 'T/m/s', 'rfRingdownTime', 30e-6, ...
                     'rfDeadTime', 100e-6, 'adcDeadTime', 20e-6);  

arg.systemGE = toppe.systemspecs('addDelays', false);   % don't add delays before creating Pulseq blocks
%arg.systemGE = toppe.systemspecs(); 

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

% system struct to be used when creating Pulseq blocks
lims = arg.system;

% Define delays to pass to 'plotseq' call
%rfdel = max(arg.system.toppe.myrfdel   = round(arg.systemSiemens.rfDeadTime*1e6);       % RF delay (us)
if 0
arg.system.toppe.daqdel     = max(round(arg.systemSiemens.adcDeadTime*1e6 - arg.system.toppe.daqdel),    0);      % ADC delay (us)
arg.system.toppe.timetrwait = max(round(arg.systemSiemens.rfRingdownTime*1e6 - arg.system.toppe.daqdel), 0);   % (us)
arg.system.toppe.myrfdel    = round(arg.systemSiemens.rfDeadTime*1e6);       % RF delay (us)
arg.system.toppe.daqdel     = round(arg.systemSiemens.adcDeadTime*1e6);      % ADC delay (us)
arg.system.toppe.timetrwait = round(arg.systemSiemens.rfRingdownTime*1e6);   % (us)
end

% initialize Pulseq sequence object
seq = mr.Sequence(lims);

% Untar files
try
	system(sprintf('tar xf %s', toppeTarFile));
catch ME
	error(ME.message);
	return;
end

% Read TOPPE scan info
geRasterTime = arg.systemGE.raster;     % GE raster time (for RF, gradients, and ADC) (sec)
max_pg_iamp  = 2^15-2;                  % max TOPPE/GE "instruction amplitude" (signed short int)
d      = toppe.utils.tryread(@toppe.readloop,           arg.loopFile);         % scanloop array
modArr = toppe.utils.tryread(@toppe.readmodulelistfile, arg.moduleListFile);   % module waveforms

% clean up TOPPE files
system('rm modules.txt scanloop.txt');
for ic = 1:length(modArr)
	system(sprintf('rm %s', modArr{ic}.fname));
end


%% Loop through scanloop.txt. Add each row as one Pulseq "block".
nt = size(d,1);    % number of startseq calls
for ii = 1:nt
	if ~mod(ii,250)
		fprintf('.');
	end

	module = modArr{d(ii,1)};

	% get waveforms and delay for one row (one startseq call)
	% rf: Gauss; gradients: Gauss/cm; tdelay: microsec
	[~, ~, ~, ~, rfwav, gxwav, gywav, gzwav, tdelay] = toppe.plotseq(ii, ii, ...
		'loopArr', d, 'mods', modArr, 'doDisplay', false, 'system', arg.systemGE);  

	% pulseq likes row vectors
	rfwav = rfwav(:)';
	gxwav = gxwav(:)';
	gywav = gywav(:)';
	gzwav = gzwav(:)';

	% convert to Pulseq units and rastertimes
	% rf:   Hz,   1us
   % grad: Hz/m, 10us
	rfwavPulseq = rf2pulseq(rfwav,geRasterTime,seq);
	gxwavPulseq = g2pulseq( gxwav,geRasterTime,seq);
	gywavPulseq = g2pulseq( gywav,geRasterTime,seq);
	gzwavPulseq = g2pulseq( gzwav,geRasterTime,seq);

	% ensure equal duration (interpolation to Pulseq rastertimes can result in unequal duration)
	% not needed?
	trf   = length(rfwavPulseq) * seq.rfRasterTime;
	tgrad = length(gxwavPulseq) * seq.gradRasterTime;
	ngradextra = ceil((trf-tgrad)/seq.gradRasterTime);
	gxwavPulseq = toppe.utils.makeevenlength( [gxwavPulseq zeros(1, ngradextra)] );
	gywavPulseq = toppe.utils.makeevenlength( [gywavPulseq zeros(1, ngradextra)] );
	gzwavPulseq = toppe.utils.makeevenlength( [gzwavPulseq zeros(1, ngradextra)] );
	tgrad  = length(gxwavPulseq) * seq.gradRasterTime;
	rfwavPulseq = [rfwavPulseq zeros(1,round((tgrad-trf)/seq.rfRasterTime))];

	% Make Pulseq gradient structs (even all zero waveforms)
	gx = mr.makeArbitraryGrad('x', gxwavPulseq, lims);
	gy = mr.makeArbitraryGrad('y', gywavPulseq, lims);
	gz = mr.makeArbitraryGrad('z', gzwavPulseq, lims);

	% bitmask indicating non-zero gradients
	hasg = 0;   
	if ~all(gxwavPulseq == 0)
		hasg = bitset(hasg,1);
	end
	if ~all(gywavPulseq == 0)
		hasg = bitset(hasg,2);
	end
	if ~all(gzwavPulseq == 0)
		hasg = bitset(hasg,3);
	end
	strArg = getArgStr(hasg);        % 'gz' or 'gx,gy,gz' or... as appropriate

	freqOffset  = d(ii,15);                         % Hz

	if module.hasRF
		phaseOffset = d(ii,12)/max_pg_iamp*pi;          % radians
		flip = pi; % d(ii,2)/max_pg_iamp*pi;

		rf = mr.makeArbitraryRf(rfwavPulseq, flip, 'FreqOffset', freqOffset, ...
			'PhaseOffset', phaseOffset, 'system', lims);

		if isempty(strArg)
			seq.addBlock(rf);
		else
			eval( sprintf( 'seq.addBlock(rf, %s)', strArg) );
		end

		if arg.debug
			clf;
			subplot(221); plot(abs(rf.signal),'r'); title(sprintf('max = %f', max(abs(rf.signal)))); ylabel('Hz');
			subplot(222); plot(angle(rf.signal),'r');  title(sprintf('max = %f', max(angle(rf.signal)))); ylabel('rad');
		end
	elseif module.hasDAQ
		phaseOffset = d(ii,13)/max_pg_iamp*pi;          % radians
		if arg.debugAdc
			% delay and shorten adc window
			nadc = 2*round(numel(gx.waveform)/4);
			delay = round(nadc/4)*seq.gradRasterTime;
			adc = mr.makeAdc(nadc, lims, 'Dwell', seq.gradRasterTime, 'delay', delay, ...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		else
			adc = mr.makeAdc(numel(gx.waveform), lims, 'Dwell', seq.gradRasterTime, 'delay', 0,...
				'freqOffset', freqOffset, 'phaseOffset', phaseOffset);
		end

		if isempty(strArg)
			seq.addBlock(adc);
		else
			eval( sprintf( 'seq.addBlock(%s, adc)', strArg) );
		end
	else
		if ~isempty(strArg)
			eval( sprintf( 'seq.addBlock(%s)', strArg) );
		end
	end

	if arg.debug
		if ~module.hasRF
			clf;
		end
		subplot(2,2,[3 4]); 
		hold on
		if bitget(hasg, 1)
			plot(gx.waveform,'r'); ylabel('Hz/m'); hold on; 
		end
		if bitget(hasg, 2)
			plot(gy.waveform,'g'); hold on;
		end
		if bitget(hasg, 3)
			plot(gz.waveform,'b'); 
		end
		if ~hasg
			clf(sfh);
		end
			
		hold off; %title(sprintf('max = %f', max([gx.waveform gy.waveform gz.waveform])));
		input('press any key to continue');
	end

	% add delay block
	if tdelay > 12;    % minimum duration of wait pulse in TOPPE
		del = mr.makeDelay(round(tdelay*1e-6,5)); %delay also needs to be in multiples of rastertimes of 10us
		seq.addBlock(del);
	end
	
end
fprintf('\n');

%seq.plot();
seq.write(arg.seqfile);

return;

	
%% convert gradient from Gauss/cm to Hz/m, and interpolate to seq.gradRasterTime
function gout = g2pulseq(g,geRasterTime,seq)
gamma = 4.2576e3;      % Hz/G
g = g * gamma * 100;   % Hz/m
T = numel(g)*geRasterTime;    % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.gradRasterTime:(T-seq.gradRasterTime);
gout = interp1(tge, g, t, 'linear', 'extrap');
return;

%% convert rf units from Gauss to Hz, and interpolate to seq.rfRasterTime
function rfout = rf2pulseq(rf,geRasterTime,seq)
gamma = 4.2576e3;       % Hz/G
rf = rf*gamma;          % Hz
T = numel(rf)*geRasterTime;   % pulse duration
tge = 0:geRasterTime:(T-geRasterTime);
t = 0:seq.rfRasterTime:(T-seq.rfRasterTime);
rfout = interp1(tge, rf, t, 'linear', 'extrap');
%L = 10; cutoff = 0.9;
%rf = interp(rf,dt/rfRasterTime,L,cutoff);      % upsample from 4us to 1us
return;

%% get gradient arguments (as string) to pass to seq.addBlock()
function argStr = getArgStr(hasg)

switch hasg
	case 0
		argStr = '';  % no gradients
	case 1
		argStr = 'gx'; 
	case 2
		argStr = 'gy'; 
	case 4
		argStr = 'gz'; 
	case 3
		argStr = 'gx, gy'; 
	case 5
		argStr = 'gx, gz'; 
	case 6
		argStr = 'gy, gz'; 
	case 7
		argStr = 'gx, gy, gz'; 
end

return;
