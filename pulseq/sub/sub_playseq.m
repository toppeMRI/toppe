function sub_playseq(modArr, loopArr, nBlocksPerTR, varargin)
% Play sequence contained in modArr and loopArr. See seq2ge.m.
%
% Example:
%  sys = toppe.systemspecs('maxSlew', 20);
%  [modArr,loopArr] = seq2ge(seq, 'system', sys, 'verbose', true);
%  sub_playseq(modArr, loopArr, 5);
%  sub_playseq(modArr, loopArr, 5, 'gradMode', 'slew');

%% parse inputs
arg.nTRskip = 0;
arg.tpause = 0;    % sec
arg.system = toppe.systemspecs();
arg.gradMode  = 'amplitude';          % 'amplitude' or 'slew'

% Substitute varargin values as appropriate
arg = toppe.utils.vararg_pair(arg, varargin);

for ii = 1:nBlocksPerTR*(1+arg.nTRskip):length(loopArr)                        
	[rf,gx,gy,gz] = sub_plotseq(modArr, loopArr, ii, ii+nBlocksPerTR-1, ...
		'nTRskip',  arg.nTRskip, ...
		'system',   arg.system, ...
		'gradMode', arg.gradMode);
	%refresh(gcf)
	pause(arg.tpause);   % helps refresh figure
end

