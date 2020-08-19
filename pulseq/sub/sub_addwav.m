%% Add waveform (column vector) to existing array. Pad with zeros as needed.
function wavs = sub_addwav(wavs, newwav)

nt   = max(size(wavs,1), size(newwav,1));
wavs = [wavs; zeros(nt-length(wavs), size(wavs,2))];
newwav = [newwav; zeros(nt-length(newwav), size(newwav,2))];
wavs = [wavs newwav];

return

