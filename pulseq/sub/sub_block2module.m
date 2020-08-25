%% Storage for groups of arbitrary gradients and rf signals
function module = sub_block2module(block, blockid, system, modnum)
%
% Initialize a module with waveforms from 'block'
%
% Inputs:
%   block            Pulseq block obtained with getBlock()
%   blockid          Pulseq block is seq.getBlock(blockid)
%   system           hardware specs, as described in ../seq2ge.m
%   modnum           Module number. Determines .mod file name.
% Output:
%   module        Struct containing module info (rf/gradient waveforms, blockEvents, etc)
%                   .duration        'duration' entry in modules.txt
%                   .hasRF           'hasRF' entry in modules.txt
%                   .hasADC          'hasADC' entry in modules.txt
%                   .rf              [nt npulses] Normalized RF waveforms 
%                   .gx              [nt npulses] Normalized Gx waveforms (same for .gy, .gz)
%                   .gy 
%                   .gz 
%                   .ofname          .mod output file name
%                   .nt              max number of samples in waveforms 
%                   .npulses         number of different sets of waveforms (e.g., size(gx.waveforms,2))
%                   .block           Pulseq block

%if length(moduleArr)+1 > system.toppe.nMaxModules
%	error(sprintf('The number of modules exceeds system.toppe.nMaxModules (%d).\nAre you sure you need that many modules?', system.toppe.nMaxModules));
%end

% Initialize with defaults
module          = struct();
module.duration = 0;
module.hasRF    = 0;
module.hasADC   = 0;
module.ofname   = sprintf('module%d.mod', modnum);
module.blockids = [];
module.modnum   = modnum;
	
% Initialize waveform arrays to []
for ax = {'rf', 'gx','gy','gz'}
	module.(ax{1}) = [];
end

module.nt = length(module.rf);   % 0
module.npulses = 0;

% Update 'module' with waveform information from 'block'.
module = sub_updatemodule(module, block, blockid, system);  

return
