function [rf, gx, gy, gz] = padwaveforms(varargin)
% function [rf, gx, gy, gz] = padwaveforms(varargin)
%
% Pad input waveform(s) with zero, and return
% arrays rf/gx/gy/gz of common size [ndat npulses],
% where ndat is multiple of 4.
%
% Input options:
%   rf            [nrf nrfpulses] Complex RF waveform, [ndat nrfpulses]
%   gx            [ngx ngxpulses]
%   gy            [ngy ngypulses]
%   gz            [ngz ngzpulses]
%
% Outputs:
%   rf            [ndat npulses], where ndat = max waveform length, and 
%                 npulses = max number of waveforms (pulses) per axis
%   gx            [ndat npulses]
%   gy            [ndat npulses]
%   gz            [ndat npulses]

import toppe.*
import toppe.utils.*

%% parse inputs
% Defaults
arg.rf = [];
arg.gx = [];
arg.gy = [];
arg.gz = [];

%arg = toppe.utils.vararg_pair(arg, varargin);
arg = vararg_pair(arg, varargin);

%% Copy input waveform to rf, gx, gy, and gz (so we don't have to carry the arg. prefix around)
fields = {'rf' 'gx' 'gy' 'gz'};
for ii = 1:length(fields)
    wavtype = fields{ii} ;   % 'rf', 'gx', 'gy', or 'gz'
    cmd = sprintf('%s = %s;', wavtype, sprintf('arg.%s', wavtype));
    eval(cmd);
end

%% Force all waveform arrays to have the same dimensions,
%% and make length multiple of 4.
ndat    = max( [size(rf,1) size(gx,1) size(gy,1) size(gz,1)] );
npulses = max( [size(rf,2) size(gx,2) size(gy,2) size(gz,2)] );
if ndat == 0
    error('At least one waveform must be non-empty');
end
if ndat > 2^15
    error(sprintf('Max waveform length is 32768 (samples) -- found %d samples', ndat));
end

% make length divisible by 4 (EPIC seems to behave best this way)
if mod(ndat, 4)
    warning('Waveform duration will be padded to 4 sample boundary.')
    ndat = ndat - mod(ndat, 4) + 4;
end

for ii = 1:length(fields)
    wavtype = fields{ii} ;                    % 'rf', 'gx', 'gy', or 'gz'
    wav = eval(fields{ii}) ;                  % [ndat npulses]

    if wavtype == 'rf' & isempty(wav)
        % Must have non-zero RF waveform to make toppe happy (even if it's not actually played out)
        wav = [zeros(1,npulses); 0.01*ones(ndat-2, npulses); zeros(1,npulses)];
    end

    % enforce equal number of rows and columns
    [nrows n2] = size(wav);

    if (nrows ~=0 & nrows < ndat) 
        warning('Padding %s with zero rows', wavtype);
    end
    if (n2 ~= 0 & n2 < npulses) 
        warning('Padding %s with zero columns', wavtype);
    end

    wav = [wav; zeros(ndat-nrows,n2)];
    wav = [wav  zeros(ndat,npulses-n2)];

    % copy to corresponding wav type (rf, gx, gy, or gz)
    cmd = sprintf('%s = %s;', wavtype, 'wav') ;
    eval (cmd);
end
    
return;

