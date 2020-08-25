%% Add waveforms from 'block' to an existing module waveform array
function module = sub_updatemodule(module, block, blockid, system)
% Add the (normalized) waveforms contained in the block numer 'blockid' 
% to the last column of module.rf/module.<grad> (if non-empty).
%
% Inputs:
%   module           see sub_block2module()
%   blockid          Pulseq block is seq.getBlock(blockid)
%   system           system struct (see above)

dt  = system.raster;   % sec

if isempty(module.blockids)
	module.blockids(1) = blockid;
else
	module.blockids(end+1) = blockid;
end

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
	module.rf = sub_addwav(module.rf, rf);         % Normalized to amplitude 1
end

% Gradients
for ax = {'gx','gy','gz'};
	ax = ax{1};
	grad = block.(ax);
	if ~isempty(grad)
		if strcmp(grad.type, 'grad')    % arbitrary shape
			% must start and end on zero
			grad.raster = grad.t(2) - grad.t(1);
			if abs(grad.waveform(1)) > 0
				grad.waveform = [0; grad.waveform];
				grad.t = [0; (grad.t + grad.raster)];
	 		end
			if abs(grad.waveform(end)) > 0
				grad.waveform = [grad.waveform; 0];
				grad.t = [grad.t; (grad.t(end) + grad.raster)];
			end

			% interpolate to GE raster time
			tge = 0:dt:(max(grad.t)-0*dt);
			wav = interp1(grad.t, grad.waveform, tge);   % interpolate to GE raster time (4us)
			wav(isnan(wav)) = 0;                         % must be due to interp1
			if wav(1) > 0
				wav = [0; wav(:)];
			end
			if wav(end) > 0
				wav = [wav(:); 0];
			end
			wav = wav/100/system.gamma;                  % Gauss/cm
		else
			% trapezoid
			wav = sub_trap2shape(grad, system.gamma);     % Gauss/cm
		end

		% check slew
		peakSlew = max(max(diff(wav,1)/(dt*1e3),[],1));
		if peakSlew > system.maxSlew
			%[peakSlew system.maxSlew blockid module.modnum size(module.(ax),2)+1]
			error(sprintf('slew rate violation at blockid %d (%s) (%.1f%%)', blockid, ax, peakSlew/system.maxSlew*100));
		end 

		% normalize and add to module 
		if norm(wav) > 0
			wav = wav/max(abs(wav(:)));
		end
		module.(ax) = sub_addwav(module.(ax), wav(:));
	end
end

% store waveform length (useful for comparing blocks)
nt = max([length(module.rf) length(module.gx) ...
	              length(module.gy) length(module.gz)]);

% pad with zeros as needed to ensure equal length of all (non-empty) waveforms
npulses = 0;
if ~isempty(module.rf)
	wav = module.rf;
	npulses = size(wav,2);
	module.rf = [wav; zeros(nt-size(wav,1), size(wav,2))];
end
gradChannels = {'gx','gy','gz'};
for ii=1:length(gradChannels)
	eval(sprintf('wav = module.%s;', gradChannels{ii}));
	if ~isempty(wav)
		npulses = size(wav,2);
		wav = [wav; zeros(nt-size(wav,1), size(wav,2))];
		eval(sprintf('module.%s= wav;', gradChannels{ii}));
	end
end

% ADC
if ~isempty(block.adc)
	module.hasADC = 1;
%	nAdc = round(block.adc.dwell/dt*block.adc.numSamples);
%	module.nt = max(nt, nAdc);
%else
%	module.nt = nt;
end

module.nt = nt;
module.npulses = npulses;

return
