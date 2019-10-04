function wf = trimzeros(wf)
% remove zeros at beginning and end of a waveform

while (wf(1) == 0)
	wf = wf(2:end);
end

while (wf(end) == 0)
	wf = wf(1:(end-1));
end
