%% Update/initialize loopStructArr 
function arg = sub_updateloopstruct(arg, block, nextblock, system, varargin)
%
% Inputs
%  arg          Either [] (empty), or a loopStruct struct (see below)
%  block        Pulseq block (getBlock()). Can be empty.
%  system       See ../seq2ge.m
%
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
	arg.rotmat = eye(3);     % 3x3 rotation matrix (added to toppev3)
end

% Substitute specified system values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

% If block is provided, fill in values as appropriate
if ~isempty(block)
	if ~isempty(block.rf)
		start_core = system.toppe.start_core_rf;    % RF module
	elseif ~isempty(block.adc) 
		start_core = system.toppe.start_core_daq;   % data acquisition module
	else
		start_core = system.toppe.start_core_grad;  % module containing only gradients (no RF or DAQ)
	end

	module.hasADC = 1;
	module.hasRF = 1;
	if ~isempty(block.delay)
		% Set textra to delay block, minus system-related overhead.
		% For an explanation of toppeDelay, see toppe_timing.pdf.
		toppeDelay = (system.toppe.timetrwait + system.toppe.timessi + start_core)*1e-6;   % sec
		siemensDelay = 0; % ?
		arg.textra = siemensDelay + block.delay.delay - toppeDelay;     % negative values will be set to zero in scanloop.txt, with a warning
		if arg.textra < 0
			keyboard
		end
		%arg.textra = max(siemensDelay + nextblock.delay.delay - toppeDelay, 0);
	end

	% If next block is a "pure" delay block (i.e., contains only a delay),
	% add that delay to textra.
	if ~isempty(nextblock) 
		if ~isempty(nextblock.delay) & isempty(nextblock.rf) & isempty(nextblock.adc) & ...
			isempty(nextblock.gx) & isempty(nextblock.gy) & isempty(nextblock.gz)
			% Set textra to delay of next block, minus system-related overhead.
         % For an explanation of toppeDelay, see toppe_timing.pdf.
         toppeDelay = (system.toppe.timetrwait + system.toppe.timessi + start_core)*1e-6;   % sec
			siemensDelay = 0; % ?
			arg.textra = arg.textra + siemensDelay + nextblock.delay.delay - toppeDelay;     % negative values will be set to zero in scanloop.txt, with a warning
			%arg.textra = max(siemensDelay + nextblock.delay.delay - toppeDelay, 0);
		end
	end

	% rf amplitude, phase, freq
	if ~isempty(block.rf)
		arg.rfamp = max(abs(block.rf.signal)/system.gamma);    % Gauss
		arg.rfphs = block.rf.phaseOffset;                      % radians
		arg.rffreq = block.rf.freqOffset;                      % Hz
	end

	% gradient amplitude
	for ax = {'gx','gy','gz'};
		ax = cell2mat(ax);
		grad = block.(ax);
		if ~isempty(grad)
			if strcmp(grad.type,'trap')
				eval( sprintf('arg.%samp = (block.%s.amplitude)/system.gamma/100;', ax, ax) );        % Gauss/cm. signed. (in moduleArr, traps are normalized and positive)
			else
				eval( sprintf('wav = block.%s.waveform/system.gamma/100;', ax) );    % Gauss/cm
				wmax = max(abs(wav));
				eval( sprintf('arg.%samp = wmax;', ax) );
				%wmax = max(wav);
				%wmin = min(wav);
				%if abs(wmin) > wmax
				%	eval( sprintf('arg.%samp = wmin') );
				%else
				%	eval( sprintf('arg.%samp = wax') );
				%end
			end
		else
		end
	end
end

return;

