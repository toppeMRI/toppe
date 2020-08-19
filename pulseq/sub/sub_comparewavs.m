%% Compare two (real-valued, positive) waveform shapes and (if same) return relative scaling
function [isSameShape] = sub_comparewavs(wav1, wav2, tol, isTrap, verbose)
% 
% Inputs
%  wav1         vector of real/complex values. If complex, only abs(wav) is considered (for now, TODO).
%  wav2         vector of real/complex values
%  tol          consider shapes to be equal if norm of difference between (normalized) waveforms < tol
%
% Outputs
%  isSameShape   boolean
%  scale         double 

if ~isreal(wav1) | ~isreal(wav2)
	wav1 = abs(wav1);
	wav2 = abs(wav2);
end

if ~exist('verbose', 'var')
	verbose = false;
end

if ~isvector(wav1) |	~isvector(wav2) | ~all(size(wav1)==size(wav2))
	error(sprintf(['input waveforms must be vectors and have the same shape (row or column)\n', ...
					'size(wav1) = (%d,%d), size(wav2) = (%d,%d)'], ...
					size(wav1,1), size(wav1,2), size(wav2,1), size(wav2,2)));
end

% normalize
wav1n = wav1/max(abs(wav1));
wav2n = wav2/max(abs(wav2));

% are shapes the same?
isSameShape = norm(wav1n-wav2n,1) < tol;          % sum(abs(difference))

% special case: trapezoids
% Pulseq produces trapezoids that are nominally scaled versions of each other, but that
% nevertheless differ enough to cause the above shape check to fail. 
% Therefore, if two trapezoids are almost identical (in shape),
% we'll assume that they can be represented as scaled versions of the one with the largest area 
% (to be determined).
if isTrap & ~isSameShape
	nslack = 10;
	if norm(wav1n-wav2n) < nslack    % allow for difference of nslack samples at peak (plateau) amplitude	
		isSameShape = true;
	end
end

if verbose
	[norm(wav1n,1) norm(wav2n,1) norm(wav1n-wav2n,1) norm(scale*wav2n,1)]
	plot([wav1n-wav2n]); title('difference after normalizing');
end

return;
