function [isValid, gmax, slewmax] = checkwaveforms(system, varargin)
% Check rf/gradient waveforms against system limits.
%
% function [isValid, gmax, slewmax] = checkwaveforms(system, varargin)
%
% Inputs:
%  system       (required) struct containing hardware specs. See systemspecs.m
%
% Options 
%  rf           rf waveform
%  gx/gy/gz     [nt npulses]. Gradient waveforms. Can be different size.
%
% Outputs
%  isValid    boolean/logical (true/false)
%  gmax       [1 3] Max gradient amplitude on the three gradient axes (x, y, z) (G/cm)
%  slewmax    [1 3] Max slew rate on the three gradient axes (x, y, z) (G/cm/ms) (G/cm/ms)

import toppe.*
import toppe.utils.*

isValid = true;

%% parse inputs

% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];

arg = toppe.utils.vararg_pair(arg, varargin);


%% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
	wavtype = fields{ii};
	cmd = sprintf('%s = %s;', wavtype, sprintf('arg.%s', wavtype));
	eval(cmd);
end

%% Zero-pad at end to equal size
[rf, gx, gy, gz] = padwaveforms('rf', rf, 'gx', gx, 'gy', gy, 'gz', gz);

%% Check against system hardware limits

axes = 'xyz';

% gradient amplitude and slew
for ii = 1:3
	eval(sprintf('g = g%s;', axes(ii))); 

    gmax(ii) = max(abs(g(:)));
    if gmax(ii) > system.maxGrad
        fprintf('Error: %s gradient amplitude exceeds system limit (%.1f%%)\n', axes(ii), gmax(ii)/system.maxGrad*100);
        isValid = false;
    end

    slewmax(ii) = 0;
    for jj = 1:size(g,2)   % loop through all waveforms (pulses)
	    slewmax(ii) = max(slewmax(ii), max(abs(diff(g(:,jj)/(system.raster*1e3)))));
    end
    if slewmax(ii) > system.maxSlew
        fprintf('Error: %s gradient slew rate exceeds system limit (%.1f%%)\n', axes(ii), slewmax(ii)/system.maxSlew*100);
        isValid = false;
    end
end

% peak rf
maxRF = max(abs(rf));
if maxRF > system.maxRF
	fprintf('Error: rf amplitude exceeds system limit (%.1f%%)\n', maxRF/system.maxRF*100);
	isValid = false;
end

%% Check PNS. Warnings if >80% of threshold.
for jj = 1:size(gx,2)   % loop through all waveforms (pulses)
    clear gtm;
    gtm(1,:) = gx(:,jj)'*1d-2;    % T/m
    gtm(2,:) = gy(:,jj)'*1d-2;
    gtm(3,:) = gz(:,jj)'*1d-2;
    [pThresh] = toppe.pns(gtm, system.gradient, 'gdt', system.raster, 'plt', false, 'print', false);
    if max(pThresh) > 80
        if max(pThresh) > 100
            warning(sprintf('PNS (%d%%) exceeds first controlled mode (100%%)!!! (waveform %d)', ...
                round(max(pThresh)), jj));
        else
            warning(sprintf('PNS(%d%%) exceeds normal mode (80%%)! (waveform %d)', ...
                round(max(pThresh)), jj));
        end
    end
end

%% do all waveforms start and end at zero?
%if any([gx(1,:) gx(end,:) gy(1,:) gy(end,:) gz(1,:) gz(end,:) rf(1,:) rf(end,:)] ~= 0)
%	fprintf('Error: all waveforms must begin and end with zero\n')
%	isValid = false;
%end

return;
