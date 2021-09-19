function [isValid, gmax, slewmax] = checkwaveforms(varargin)
% Check rf/gradient waveforms against system limits.
%
% function isValid = checkwaveforms(varargin)
%
% Inputs:
%  system       (required) struct containing hardware specs. See systemspecs.m
%
% Options 
%  rf           rf waveform
%  gx/gy/gz     gradient waveform
%  rfUnit       Gauss (default) or mT
%  gradUnit     Gauss/cm (default) or mT/m
%
% Outputs
%  isValid    boolean/logical (true/false)
%  gmax       [1 3] Max gradient amplitude on the three gradient axes (x, y, z) (G/cm)
%  slewmax    [1 3] Max slew rate on the three gradient axes (x, y, z) (G/cm/ms) (G/cm/ms)

import toppe.*
import toppe.utils.*

%% parse inputs

% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];
arg.rfUnit   = 'Gauss';
arg.gradUnit = 'Gauss/cm';
arg.system = [];

arg = toppe.utils.vararg_pair(arg, varargin);

if isempty(arg.system)
    error('Missing system argument');
end

system = arg.system;

%% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
	wavtype = fields{ii};
	cmd = sprintf('%s = %s;', wavtype, sprintf('arg.%s', wavtype));
	eval(cmd);
end

%% Convert input waveforms and system limits to Gauss and Gauss/cm
if strcmp(arg.rfUnit, 'mT')
	rf = rf/100;   % Gauss
end
if strcmp(arg.gradUnit, 'mT/m')
	gx = gx/10;    % Gauss/cm
	gy = gy/10;
	gz = gz/10;
end

if strcmp(system.rfUnit, 'mT')
	system.maxRf = system.maxRf/100;      % Gauss
end
if strcmp(system.gradUnit, 'mT/m')
	system.maxGrad = system.maxGrad/10;   % Gauss/cm
end
if strcmp(system.slewUnit, 'T/m/s')
	system.maxSlew = system.maxSlew/10;   % Gauss/cm/msec
end

%% Check against system hardware limits
isValid = true;

grads = 'xyz';

% gradient amplitude and slew
for ii = 1:3
	eval(sprintf('g = g%s;', grads(ii))); 
    if isempty(g)
        gmax(ii) = 0;
        slewmax(ii) = 0;
        continue; 
    end

    gmax(ii) = max(abs(g(:)));
    if gmax(ii) > system.maxGrad
        fprintf('Error: %s gradient amplitude exceeds system limit (%.1f%%)\n', grads(ii), gmax(ii)/system.maxGrad*100);
        isValid = false;
    end

    slewmax(ii) = 0;
    for jj = 1:size(g,2)   % loop through all waveforms (pulses)
	    slewmax(ii) = max(slewmax(ii), max(abs(diff(g(:,jj)/(system.raster*1e3)))));
    end
    if slewmax(ii) > system.maxSlew
        fprintf('Error: %s gradient slew rate exceeds system limit (%.1f%%)\n', grads(ii), slewmax(ii)/system.maxSlew*100);
        isValid = false;
    end
end

% peak rf
maxRf = max(abs(rf));
if maxRf > system.maxRf
	fprintf('Error: rf amplitude exceeds system limit (%.1f%%)\n', maxRf/system.maxRf*100);
	isValid = false;
end

%% Check PNS. Warnings if >80% of threshold.
for ii = 1:3
	eval(sprintf('g = g%s;', grads(ii))); 
    if isempty(g)
        continue; 
    end
    nwavs = size(g,2);
    for jj = 1:size(g,2)   % loop through all waveforms (pulses)
        clear gtm;
        gtm(1,:) = g(:,jj)'*1d-2;    % T/m
        gtm(2,:) = g(:,jj)'*1d-2;
        gtm(3,:) = g(:,jj)'*1d-2;
        [pThresh] = toppe.pns(gtm, system.gradient, 'gdt', system.raster, 'plt', false, 'print', false);
        if max(pThresh) > 80
            if max(pThresh) > 100
                warning(sprintf('Slew (%d%%) exceeds first controlled mode (100%%)!!! (g%s, waveform %d)', ...
                    round(max(pThresh)), grads(ii), jj));
            else
                warning(sprintf('Slew (%d%%) exceeds normal mode (80%%)! (g%s, waveform %d)', ...
                    round(max(pThresh)), grads(ii), jj));
            end
        end
	end
end

%% Is (max) waveform duration on a 4 sample (16us) boundary?
ndat = max( [size(rf,1) size(gx,1) size(gy,1) size(gz,1)] );
if mod(ndat, 4)
	fprintf('Error: waveform duration must be on a 4 sample (16 us) boundary.');
	isValid = false;
end

%% do all waveforms start and end at zero?
for ii = 1:3
	eval(sprintf('if isempty(g%s); g%s = 0; end', grads(ii), grads(ii)));
end
if isempty(rf)
	rf = 0;
end
if any([gx(1,:) gx(end,:) gy(1,:) gy(end,:) gz(1,:) gz(end,:) rf(1,:) rf(end,:)] ~= 0)
	fprintf('Error: all waveforms must begin and end with zero\n')
	isValid = false;
end

return;
