function [moduleArr loopStructArr] = seq2ge(seqarg, varargin)
% function seq2ge(seqarg, varargin)
%
% Pulseq to TOPPE file conversion.
%
% This script creates a .tar file containing the following:
%   *.mod:            One .modfile corresponds to one "unique" block (see below)
%   modules.txt       List of .modfiles, and flags indicating whether the .wav file is RF/ADC/(gradients only)
%   scanloop.txt      Sequence of instructions for the entire scan (waveform amplitudes, ADC instructions, etc)
%
% Inputs:
%   seqarg            Either a Pulseq file name, or an mr.Sequence object.
% Input options:
%  tarfile            Output file name. Default: 'out.tar'
%  system             struct containing GE and TOPPE system specs. See +toppe/systemspecs.m.
%  verbose            true (default) or false
%
% Output:
%  TOPPE tar file containing scanloop.txt, modules.txt, and .mod files
%
% Usage examples:
%   >> seq2ge('myseqfile.seq');
%
%   >> lims = toppe.systemspecs('maxSlew', 130, 'slewUnit', 'T/m/s');
%   >> seq2ge('myseqfile.seq', 'system', lims);
%
%   >> seq = mr.Sequence();
%   >> seq.read('myseqfile.seq');
%   >> seq2ge(seq);
%
% See https://toppemri.github.io/ for more info on TOPPE.
%
% $Id: seq2ge.m,v 1.26 2018/11/05 18:15:57 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/pulseq/seq2ge.m,v $

%% parse inputs
% Defaults
arg.tarfile = 'out.tar';
arg.system  = toppe.systemspecs();
arg.verbose = false;

%  systemSiemens      struct containing Siemens system specs. 
%                        .rfRingdownTime     Default: 30e-6   (sec)
%                        .rfDeadTime         Default: 100e-6  (sec)
%                        .adcDeadTime        Default: 20e-6   (sec)
%arg.systemSiemens = struct('rfRingdownTime', 30e-6, 'rfDeadTime', 100e-6,'adcDeadTime', 20e-6);

% Substitute specified system values as appropriate (from MIRT toolbox)
arg = toppe.utils.vararg_pair(arg, varargin);

switch arg.system.toppe.version
	case 'v2' 
		nCols = 16;   % number of columns in scanloop.txt
	otherwise
		error('Please use TOPPE v2');
end

%% Get seq object
if isa(seqarg, 'char')
	seq = mr.Sequence();
	seq.read(seqarg);
else
	if ~isa(seqarg, 'mr.Sequence')
		error('Input not an mr.Sequence object');
	end
	seq = seqarg;
end

%% Loop through blocks and build 'moduleArr' and 'loopStructArr'

% 'moduleArr' struct array
% Find blocks that are unique in terms of waveforms and timing (i.e., waveform amplitudes, RF/ADC phase, etc can differ),
% and fills 'moduleArr' struct array accordingly. 
% Each entry of 'moduleArr' is a struct containing all waveforms belonging to one module (.mod file), and other module info.
% The usage of the word "module" here is consistent with its use in TOPPE.

% 'loopStructArr' struct array
% Each entry contains information needed to fill out one row of scanloop.txt.

if arg.verbose
	fprintf('Filling moduleArr struct array, and loopStructArr array.\n' );
end

blockEvents = cell2mat(seq.blockEvents);
blockEvents = reshape(blockEvents, [6, length(seq.blockEvents)]).'; %hardcoded for as long as Pulseq does not include another element.

% First entry in 'moduleArr' struct array
ib = 1;
block = seq.getBlock(ib);
if ~isempty(block.delay)
	error('First block can''t contain a delay. Edit the .seq file.');
end
moduleArr(ib) = addModule([], block, arg.system);

% First entry in 'loopStructArr' struct array (first block is by definition a module)
nextblock = seq.getBlock(ib+1);   % needed to set 'textra' in scanloop.txt
loopStructArr(ib) = updateLoopStruct([], block, nextblock, arg.system, 'mod', ib);

% data frames (in Pfile) are stored using indeces 'slice', 'echo', and 'view' 
sl = 1;
view = 1;
echo = 0; 
adcCount = 0;

for ib = 2:length(seq.blockEvents)
	if ~mod(ib, 1000)
		fprintf('.');
	end

	block = seq.getBlock(ib);

	if ib < length(seq.blockEvents)
		nextblock = seq.getBlock(ib+1);  % used to set textra column in scanloop.txt
	else
		nextblock = [];
	end

	if ~isempty(block.delay)
		continue;  % Done, move on to next block (delays are accounted for in 'textra' in the previous row in scanloop.txt)
	end

	% create a TOPPE module struct from current Pulseq block
	modCandidate = addModule([], block, arg.system);

	% set slice/echo/view indeces
	% view = 1, ..., system.maxView
	% sl   = 1, ..., system.maxSlice
	if modCandidate.hasADC
		view = mod(adcCount, arg.system.maxView) + 1;
		sl   = floor(adcCount/arg.system.maxView) + 1;
		if sl > arg.system.maxSlice;
			error(sprintf('max number of slices ecxeeded (%d)', arg.system.maxSlice));
		end
		echo = floor(adcCount/(arg.system.maxView*arg.system.maxSlice));
		if echo > arg.system.maxEcho
			error('Exceeded system.maxEcho -- acquisition has too many acquired frames');
		end
		%fprintf('ib: %d, view: %d, sl: %d, echo: %d\n', ib, view, sl, echo);

		adcCount = adcCount+1;
	end

	% Does one of the existing modules (elements of moduleArr) have the same length waveform, 
	% the same non-empty rf/gx/gy/gz, and the same value of 'hasADC', as modCandidate?
	isUnique = 1; 
	for ic = 1:length(moduleArr)
		if (moduleArr(ic).nt == modCandidate.nt ...
			& isempty(moduleArr(ic).rf.waveforms) == isempty(modCandidate.rf.waveforms) ...
			& isempty(moduleArr(ic).gx.waveforms) == isempty(modCandidate.gx.waveforms) ...
			& isempty(moduleArr(ic).gy.waveforms) == isempty(modCandidate.gy.waveforms) ...
			& isempty(moduleArr(ic).gz.waveforms) == isempty(modCandidate.gz.waveforms) ...
			& moduleArr(ic).hasRF  == modCandidate.hasRF ...
			& moduleArr(ic).hasADC == modCandidate.hasADC ...
			)
			isUnique = 0;
			break;   % break out of 'for ic' loop. 'ic' now has the value of a module we'll reuse
		end
	end

	if isUnique
		% We found a unique block, so add it as a new module
		if arg.verbose
			fprintf('\tFound new module at block %d\n', ib);
		end
		moduleArr = addModule(moduleArr, block, arg.system);
		loopStructArr(ib) = updateLoopStruct([], block, nextblock, arg.system, ...
			'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', length(moduleArr));
		continue; % done, so move on to next block
	end

	% modCandidate is not unique but has the same waveform length as an existing module, 
   % so now check to see if we can find an existing set of waveforms (a column in moduleArr(ic).rf/gx/gy/gz)
	% with the same shape so we can 'reuse' it by simply scaling.
	% Typically, this will be an rf waveform or trapezoidal phase-encode gradient.
	% Note: If modCandidate contains more than one waveform type (rf/gx/gy/gz) and a match is not found
	% for ALL waveform types, then we can't reuse waveform.

	canReuse = 0;

	tol = 1e-10;  % Shape is deemed equal if sum(abs(difference)) < tol

	% rf
	%if ~isempty(block.rf) & moduleArr(ic).hasRF
	if modCandidate.hasRF
		for iwav = 1:size(moduleArr(ic).rf.waveforms, 2)
			try
				err = norm( moduleArr(ic).rf.waveforms(:,iwav) - modCandidate.rf.waveforms, 1);
			catch
				keyboard
			end
			if norm( moduleArr(ic).rf.waveforms(:,iwav) - modCandidate.rf.waveforms, 1) < tol
				% Found an rf waveform with the same shape
				canReuse = 1;
				iWavReuse = iwav;
				break;
			end
		end
	end

	% gradients
	if ~canReuse
		gradChannels={'gx','gy','gz'};
		for j=1:length(gradChannels)
			ax = gradChannels{j};
			grad = block.(ax);
			if ~isempty(grad)
				eval(sprintf('npulses = size(moduleArr(ic).%s.waveforms, 2);', ax));
				for iwav = 1:npulses
					eval(sprintf('wav1 = moduleArr(ic).%s.waveforms(:,iwav);', ax));
					eval(sprintf('wav2 = modCandidate.%s.waveforms;', ax));
					wavdiff = norm(wav1-wav2, 1); 
					if wavdiff < tol
						% Found a gradient waveform (probably a trapezoid) with the same shape.
						canReuse = 1;
						iWavReuse = iwav;
						break;
					end
				end
			end
		end
	end

	if canReuse
		% We found a set of RF/gradient waveforms we can reuse, 
		% so set 'mod' and 'wavnum' (waveform array column index) accordingly.
		if isempty(iWavReuse)
			error(sprintf('Error: empty iwav at ib=%d', ib));
		end
		loopStructArr(ib) = updateLoopStruct([], block, nextblock, arg.system, ...
			'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic, 'wavnum', iWavReuse);
	else
		% can't reuse an existing waveform, so add this waveform to moduleArr(ic)
		moduleArr(ic) = updateModule(moduleArr(ic), block, arg.system);
		loopStructArr(ib) = updateLoopStruct([], block, nextblock, arg.system, ...  %'mod', ic);
			'dabmode', 1, 'slice', sl, 'echo', echo, 'view', view, 'mod', ic);
	end

end

if arg.verbose
	fprintf(' done\n');
end


%% Write each module to a .mod file, and create modules.txt

if arg.verbose
	fprintf(1, 'Writing .mod files and modules.txt... ');
end

% write modules.txt header
fid = fopen('modules.txt','w');
fprintf(fid,'Total number of unique cores\n');
fprintf(fid,'%d\n', length(moduleArr));
fprintf(fid,'wavfile_name    duration (us)     has_RF?     has_ADC?\n');

% loop through moduleArr
for ic = 1:length(moduleArr)

	hasadc = moduleArr(ic).hasADC;
	hasrf  = moduleArr(ic).hasRF;

	if hasrf & hasadc
		error('Cannot transmit RF and acquire data in same block. Redesign the .seq file.');
	end

	% Find max rf amplitude in module and write the module using that value.
	% Then, if instruction amp = 32766 in scanloop.txt, the waveform should play out at the desired max rf.

	rfmax = 0;
	if hasrf
		for ib = 1:length(loopStructArr)
			if loopStructArr(ib).mod == ic
				rfmax = max(loopStructArr(ib).rfamp, rfmax);
			end
		end
	end

	% find max gradient amplitudes
	gxmax = 0;
	gymax = 0;
	gzmax = 0;
	for ib = 1:length(loopStructArr)
		if loopStructArr(ib).mod == ic
			gxmax = max(loopStructArr(ib).gxamp, gxmax);
			gymax = max(loopStructArr(ib).gyamp, gymax);
			gzmax = max(loopStructArr(ib).gzamp, gzmax);
		end
	end

	%fprintf('module %d: rfmax: %.3f, gxmax:%.2f, gymax:%.2f, gzmax:%.2f\n', ic, rfmax, gxmax, gymax, gzmax);

	% waveforms (ok if empty)
	rf = moduleArr(ic).rf.waveforms * rfmax;
	gx = moduleArr(ic).gx.waveforms;
	gy = moduleArr(ic).gy.waveforms;
	gz = moduleArr(ic).gz.waveforms;

	if arg.verbose
		fprintf('Creating .mod file number %d... ', ic);
	end

	% force rf waveform to be non-empty to avoid error in writemod (does no harm; will not be played out)
	if isempty(rf) 
		rf = 0.001*ones(moduleArr(ic).nt, 1);
	end

	try
		warning('off');    % don't show message about padding waveforms
		toppe.writemod('system', arg.system, 'rf', rf, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', moduleArr(ic).ofname); 
		warning('on');
	catch ME
		error(sprintf('Error in writemod:\n%s', ME.message));
	end
		
	if arg.verbose
		fprintf('success\n');
		figure; toppe.plotmod(moduleArr(ic).ofname); subplot(331); title(sprintf('module %d', ic));
	end

	% update entry in modules.txt
	fprintf(fid,'%s\t%d\t%d\t%d\n', moduleArr(ic).ofname, 0, hasrf, hasadc);	
end
fclose(fid);

if arg.verbose
	fprintf('done. Created %d .mod files.\n', ic);
end


%% Write scanloop.txt, which specifices the scan sequence (along with modules.txt and the .mod files).

% load .mod files
mods = toppe.utils.tryread(@toppe.readmodulelistfile, 'modules.txt');

if arg.verbose
	fprintf('Writing scanloop.txt... ');
end

max_pg_iamp = 2^15-2;  
ia_th = max_pg_iamp;

d = zeros(length(loopStructArr), nCols);    % will be pruned later

% loop through blocks
sendTextraWarning = true;
il = 0;   % index into loopStructArr struct (not the same as 'ib' next b/c of delay blocks
for ib = 1:length(loopStructArr)

	if isempty(loopStructArr(ib).mod)
		% skip delay blocks
		continue;
	end

	il = il + 1;

	% module and wave number
	ic   = loopStructArr(ib).mod;        % module to use for this block
	iwav = loopStructArr(ib).wavnum;     % points to RF/gradient waveform to be used [TODO: rfwavnum, gxwavnum]
	d(il, 1) = ic;
try
	d(il, 16) = iwav;
catch
	error(sprintf('error at ib=%d, il=%d', ib, il));
end

	% rf 'instruction' amplitude (even int16)
	%if mods{ic}.hasRF
	if mods{ic}.hasRF
		rf = mods{ic}.rf;                % Gauss
		peak = max(abs(rf(:)));          % Gauss
		if peak > 0
			d(il, 2) = 2*floor(max_pg_iamp*loopStructArr(ib).rfamp/peak/2);
			2*round(max_pg_iamp*loopStructArr(ib).rfamp/peak/2);
		end
	end

	% theta waveform instruction amplitude
	d(il, 3) = max_pg_iamp;      

	% gradient instruction amplitudes
	if ~isempty(moduleArr(ic).gx.waveforms)
		peak = max(abs(moduleArr(ic).gx.waveforms(:)));
		if peak > 0
			d(il, 4) = 2*floor(max_pg_iamp * loopStructArr(ib).gxamp / peak /2);
		end
	end
	if ~isempty(moduleArr(ic).gy.waveforms)
		peak = max(abs(moduleArr(ic).gy.waveforms(:)));
		if peak > 0
			d(il, 5) = 2*floor(max_pg_iamp * loopStructArr(ib).gyamp / peak /2);
		end
	end
	if ~isempty(moduleArr(ic).gz.waveforms)
		peak = max(abs(moduleArr(ic).gz.waveforms(:)));
		if peak > 0
			d(il, 6) = 2*floor(max_pg_iamp * loopStructArr(ib).gzamp / peak /2);
		end
	end

	% slice, echo, view
	d(il, 7) = loopStructArr(ib).slice;
	d(il, 8) = loopStructArr(ib).echo;
	d(il, 9) = loopStructArr(ib).view;

	% dabmod
	d(il, 10) = loopStructArr(ib).dabmode;

	% in-plane (kx, ky) rotation (+/-max_pg_iamp = +/-pi radians)
	d(il, 11) = 2*floor(max_pg_iamp*loopStructArr(ib).phi/pi/2);;

	% RF transmit and receive phase (+/-max_pg_iamp = +/-pi radians)
	d(il, 12) = 2*floor(max_pg_iamp*loopStructArr(ib).rfphs/pi/2);;
	d(il, 13) = 2*floor(max_pg_iamp*loopStructArr(ib).recphs/pi/2);;

	% textra 
	if loopStructArr(ib).textra < 0
		sendTextraWarning = true;
	end
	textra = max(loopStructArr(ib).textra, 0);            % sec
	d(il, 14) = 2*round(textra*1e6/2);                    % microseconds

	% rf transmit frequency (Hz, integer)
	d(il, 15) = 2*round(loopStructArr(ib).rffreq/2);
end

% write to scanloop.txt
toppe.loop2txt(d(1:il,:));

if sendTextraWarning
	fprintf('\nWarning: .seq sequence timing too tight -- ''textra'' column set to zero in one or more scanloop.txt entries.');
end

if arg.verbose
	fprintf(' done\n');
end


%% Put TOPPE files in a .tar file (for convenience)
system(sprintf('tar cf %s modules.txt scanloop.txt', arg.tarfile));
for ic = 1:length(moduleArr)
	system(sprintf('tar rf %s %s', arg.tarfile, moduleArr(ic).ofname));
end

% list archive contents
if arg.verbose
	fprintf('\nCreated %s containing the following files:\n', arg.tarfile);
	system(sprintf('tar tf %s', arg.tarfile));
end

% clean up
system('rm modules.txt scanloop.txt');
for ic = 1:length(moduleArr)
	system(sprintf('rm %s', moduleArr(ic).ofname));
end

if ~arg.verbose
	fprintf('\n');
end

fprintf('Remember to rename one of the .mod files to ''tipdown.mod'', and another to ''readout.mod''\n');

return

%% End of main script

		
%% Storage for groups of arbitrary gradients and rf signals
function moduleArr = addModule(moduleArr, block, system)
%
% Inputs:
%   moduleArr     Empty ([]), or a struct array containing module info (rf/gradient waveforms, blockEvents, etc)
%                 Fields:
%                   .duration        'duration' entry in modules.txt
%                   .hasRF           'hasRF' entry in modules.txt
%                   .hasADC          'hasADC' entry in modules.txt
%                   .rf.waveforms    [nt npulses] Normalized RF waveforms 
%                   .gx.waveforms    [nt npulses] Normalized Gx waveforms (same for .gy, .gz)
%                   .gy.waveforms    [nt npulses] Normalized Gy waveforms (same for .gy, .gz)
%                   .gz.waveforms    [nt npulses] Normalized Gz waveforms (same for .gy, .gz)
%                   .ofname          .mod output file name
%                   .nt              max number of samples in waveforms 
%   block            Pulseq block obtained with getBlock()
%   system           hardware specs, as described above

if length(moduleArr)+1 > system.toppe.nMaxModules
	error(sprintf('The number of modules exceeds system.toppe.nMaxModules (%d).\nAre you sure you need that many modules?', system.toppe.nMaxModules));
end

% Initialize with defaults
module          = struct();
module.duration = 0;
module.hasRF    = 0;
module.hasADC   = 0;
module.ofname   = sprintf('module%d.mod', length(moduleArr)+1);
	
% Initialize waveform arrays to []
module.rf.waveforms = []; 
gradChannels={'gx','gy','gz'};
for j=1:length(gradChannels)
	ax = gradChannels{j};        
	module.(ax).waveforms = [];
end
module.nt = length(module.rf.waveforms);   % 0

% Update 'module' with waveform information from 'block'.
module = updateModule(module, block, system);  

if isempty(moduleArr)
	% initialize 'moduleArr' struct array
	moduleArr = module;
else
	% add to existing struct array
	moduleArr(end+1) = module;
end

return


%% Add waveforms from 'block' to module waveform array
function module = updateModule(module, block, system)
% Add the waveforms contained in 'block' to the last column of
% module.rf.waveforms/module.<grad>.waveforms (if non-empty).
% Also set module.hasADC and module.hasRF.
%
% Inputs:
%   module           One member of 'moduleArr' struct array.
%   block            Pulseq block obtained with getBlock()
%   system           system struct (see above)

dt  = system.raster;   % sec

% RF
if ~isempty(block.rf)
	module.hasRF = 1;

	% interpolate to GE raster time (4us)
	rf = downsample(block.rf.signal, round(dt/1e-6));     % downsample from 1us to 4us (GE raster time)
	%if ( length(block.rf.signal) > 24 )
	%else
	%	rf = block.rf.signal;                               % assumed to be already decimated in terms of dt_ge
	%	warning('rf waveform is < 24 points and will not be interpolated to GE raster time');
	%end

	% Normalize and add to waveforms (loopStructArr array will contain rf amplitude for each block)
	rf = rf/max(abs(rf(:)));
	module.rf.waveforms = [module.rf.waveforms rf];         % Normalized to amplitude 1
end

% Gradients
gradChannels = {'gx','gy','gz'};
for j=1:length(gradChannels)
	ax = gradChannels{j};        
	grad = block.(ax);
	if ~isempty(grad)
		if strcmp(grad.type, 'grad')
			% arbitrary gradient
			tge = 0:dt:(max(grad.t)-dt);
			wav = interp1(grad.t, grad.waveform, tge);   % interpolate to GE raster time (4us)
			wav(isnan(wav)) = 0;                         % must be due to interp1
			wav = wav/100/system.gamma;                  % Gauss/cm
			maxSlew = max(abs(diff(wav)/4-3));
			if maxSlew > system.maxSlew
%				error(sprintf('A %s arbitrary gradient exceeds slew rate limit', ax));
			end	
			module.(ax).waveforms = [module.(ax).waveforms wav(:)];   % ok to have non-even number of samples here; writemod.m will fix that.
		else
			% trapezoid
			gamp = grad.amplitude/100/system.gamma;                         % Gauss/cm
			dgmax = system.maxSlew*dt;                                      % max gradient change per waveform sample (G/cm)
			newwav = [ linspace(0, gamp, ceil(grad.riseTime/dt))  ...
						gamp*ones(1, floor(grad.flatTime/dt)) ... 
						linspace(gamp, 0, ceil(grad.fallTime/dt)) ]';
			maxSlew = max(abs(diff(newwav)/4e-3));
			if maxSlew > system.maxSlew
%				error(sprintf('A %s trapezoid gradient exceeds slew rate limit (%.1f%%)', ax, maxSlew/system.maxSlew*100));
			end	
			% ensure that length is consistent with existing waveform array
			wavs = module.(ax).waveforms;
			nt   = max(length(wavs), length(newwav));   % new waveform length
			wavs = [wavs; zeros(nt-length(wavs), size(wavs,2))];
			newwav = [newwav; zeros(nt-length(newwav), 1)];

			% add new waveform to array
			module.(ax).waveforms = [wavs newwav];
		end
	end
end

% store waveform length (useful for comparing blocks)
nt = max([length(module.rf.waveforms) length(module.gx.waveforms) ...
	              length(module.gy.waveforms) length(module.gz.waveforms)]);

% pad waveforms with zeros as needed to achieve equal length
%if ~isempty(module.rf.waveforms)
%	wav = module.rf.waveforms;
%	module.rf.waveforms = [wav; zeros(nt-size(wav,1), size(wav,2))];
%end
%gradChannels = {'gx','gy','gz'};
%for ii=1:length(gradChannels)
%	eval(sprintf('wav = module.%s.waveforms;', gradChannels{ii}));
%	if ~isempty(wav)
%		wav = [wav; zeros(nt-size(wav,1), size(wav,2))];
%		eval(sprintf('module.%s.waveforms = wav;', gradChannels{ii}));
%	end
%end

% ADC
if ~isempty(block.adc)
	module.hasADC = 1;
	nAdc = round(block.adc.dwell/dt*block.adc.numSamples);
	module.nt = max(nt, nAdc);
else
	module.nt = nt;
end

return


%% Update/initialize loopStructArr 
function arg = updateLoopStruct(arg, block, nextblock, system, varargin)
%
% Inputs
%  arg          Either [] (empty), or a loopStruct struct (see below)
%  block        Pulseq block (getBlock()). Can be empty.
%  system       See above
% Options
%  see code below

% Initialize loopStruct struct
if isempty(arg)
	% Defaults
	arg.mod   = 1;           % module number. Positive integer (starts at 1).
	arg.rfamp = 0;           % rf waveform (rho) amplitude (Gauss)
	arg.gxamp = 0;           % Gx waveform amplitude       (Gauss/cm)
	arg.gyamp = 0;           % Gy waveform amplitude
	arg.gzamp = 0;           % Gz waveform amplitude
	arg.slice = 1;           % data storage ’slice’ index. Positive integer (starts at 1).
	arg.echo  = 0;           % data storage ’echo’ index.  Non-negative integer (starts at 0).
	arg.view  = 1;           % data storage ’view’ index.  Positive integer (starts at 1).
	arg.dabmode = 1;         % turn on/off data acquisition (1/0)
	arg.phi    = 0;          % in-plane (x-y) rotation angle (radians)
	arg.rfphs  = 0;          % RF transmit phase (radians)
	arg.recphs = 0;          % receive phase (radians)
	arg.textra = 0;          % time added to end of module (sec)
	arg.rffreq = 0;          % RF transmit frequency offset (Hz)
	arg.wavnum = 1;          % waveform number (rf/grad waveform array column index). Non-zero positive integer.
end

% Substitute specified system values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

% If block is provided, fill in values as appropriate
if ~isempty(block)
	if ~isempty(nextblock) 
		if ~isempty(nextblock.delay)
			% Set textra to delay of next block, minus system-related overhead.
         % For an explanation of toppeDelay, see toppe_timing.pdf.
         toppeDelay = (system.toppe.timetrwait + system.toppe.timessi + system.toppe.start_core)*1e-6;   % sec
			siemensDelay = 0; % ?
			arg.textra = siemensDelay + nextblock.delay.delay - toppeDelay;     % negative values will be set to zero in scanloop.txt, with a warning
			%arg.textra = max(siemensDelay + nextblock.delay.delay - toppeDelay, 0);
		end
	end

	if ~isempty(block.rf)
		arg.rfamp = max(abs(block.rf.signal)/system.gamma);    % Gauss
		arg.rfphs = block.rf.phaseOffset;                      % radians
		arg.rffreq = block.rf.freqOffset;                      % Hz
	end

	gradChannels={'gx','gy','gz'};
	for j=1:length(gradChannels)
		ax = gradChannels{j};
		grad = block.(ax);
		if ~isempty(grad)
			if strcmp(grad.type,'trap')
				eval( sprintf('arg.%samp = abs(block.%s.amplitude)/system.gamma/100;', ax, ax) );        % Gauss/cm
			else
				eval( sprintf('arg.%samp = max(abs(block.%s.waveform))/system.gamma/100;', ax, ax) );    % Gauss/cm
			end
		else
		end
	end
end

return;



