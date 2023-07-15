function [rf, gx, gy, gz] = padwaveforms(varargin)
% function [rf, gx, gy, gz] = padwaveforms(varargin)
%
% Return arrays rf/gx/gy/gz of common size [ndat npulses],
% where ndat = max waveform length, and 
% npulses = max number of waveforms (pulses) per channel/axis
%
% Input options:
%   rf            [nrf nrfpulses] Complex RF waveform
%   gx            [ngx ngxpulses]
%   gy            [ngy ngypulses]
%   gz            [ngz ngzpulses]
%
% Outputs:
%   rf            [ndat npulses], where ndat = max([nrf ngx ngy ngz]) and
%                 npulses = max([nrfpulses ngxpulses ngypulses ngzpulses])
%   gx            [ndat npulses]
%   gy            [ndat npulses]
%   gz            [ndat npulses]

import toppe.*
import toppe.utils.*

% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];

% replace defaults with user inputs
%arg = toppe.utils.vararg_pair(arg, varargin);
arg = vararg_pair(arg, varargin);

% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
    channel = fields{ii} ;   % 'rf', 'gx', 'gy', or 'gz'
    cmd = sprintf('%s = %s;', channel, sprintf('arg.%s', channel));
    eval(cmd);
end

ndat    = max( [size(rf,1) size(gx,1) size(gy,1) size(gz,1)] );
npulses = max( [size(rf,2) size(gx,2) size(gy,2) size(gz,2)] );
if ndat == 0
    error('At least one waveform must be non-empty');
end
if ndat > 2^15
    warning(sprintf('waveform length is %d samples', ndat));
end

for ii = 1:length(fields)
    channel = fields{ii} ;                    % 'rf', 'gx', 'gy', or 'gz'
    wav = eval(fields{ii}) ;                  % [ndat npulses]

    if isempty(wav)
        wav = zeros(ndat, npulses);
    end

    if channel == 'rf'
        if isempty(wav) | norm(wav) == 0
            % Must have non-zero RF waveform to make toppe happy (even if it's not actually played out)
            wav = [zeros(1, npulses); 0.01*ones(ndat-2, npulses); zeros(1, npulses)];
        end
    end

    [nrows ncols] = size(wav);
    wavnew = zeros(ndat, npulses);

    for ic = 1:ncols
        wavnew(:,ic) = [wav(:,ic); wav(end,ic)*ones(ndat-nrows, 1)];
    end

    cmd = sprintf('%s = wavnew;', channel) ;
    eval (cmd);
end
    
return;

