function sub_playseq(modArr, loopArr, nBlocksPerTR, nTRskip, tpause)
% Play sequence contained in modArr and loopArr. See seq2ge.m.
% Example
%  [modArr,loopArr] = seq2ge(seq, 'system', lims, 'verbose', false);
%  sub_playseq(modArr, loopArr, 5, 4);

if ~exist('nTRskip', 'var')
	nTRskip = 0;
end
if ~exist('tpause', 'var')
	tpause = 0.05;   % sec
end

for ii = 1:nBlocksPerTR*(1+nTRskip):length(loopArr)                        
	[rf,gx,gy,gz] = sub_plotseq(modArr,loopArr,ii,ii+nBlocksPerTR-1);
	%refresh(gcf)
	pause(tpause);   % helps refresh figure
end

