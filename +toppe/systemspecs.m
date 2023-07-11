function sys = systemspecs(varargin)
% Create struct containing scanner hardware specs and related info.
%
% function sys = systemspecs(varargin)
% 
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
%
% Times are in microseconds, except raster (sec)


%% Defaults

% Global constants
sys.raster = 4e-6;              % sec. Raster time for gradient and RF waveforms.
sys.gamma  = 42.576e6;          % Hz/T

% Scanner-specific settings
sys.B0 = 3.0;                     % field strength (T)
sys.gradient = 'xrm';             % gradient coil
sys.psd_grd_wait = 100;           % ADC/gradient delay (us).
sys.psd_rf_wait = 100;            % rf/gradient delay (us)
sys.segmentDeadTime = 12;         % Dead time before start of block group, equals RUP_GRD(9us)
sys.segmentRingdownTime = 104;    % Delay at end of block group, equals 4us + timssi. 
sys.forbiddenEspRange = [0.41 0.51];    % (ms) Forbidden echo spacings (mechanical resonance). See /usr/g/bin/epiesp.dat
sys.tminwait   = 12;              % minimum duration of wait pulse in EPIC code (us)

% Design choices (need not equal scanner limits)
sys.maxGrad  = 4;      % Gauss/cm
sys.maxSlew  = 15;     % Gauss/cm/ms
sys.maxRF    = 0.15;   % Not sure what the hardware limit is here
sys.rfDeadTime = 72;         % us. Must be >= 72us
sys.rfRingdownTime = 54;     % us. Must be >= 54us
sys.adcDeadTime = 40;        % us. Must be >= 40us

% The following determine the slice/echo/view indexing in the data file
sys.maxSlice = 2048;           % max dabslice. UI won't allow more than this to be entered
sys.maxView  = 600;            % not sure what limit is here
sys.maxEcho  = 16;             % determined empirically


%% Substitute specified system values
sys = toppe.utils.vararg_pair(sys, varargin);

sys.gradient = lower(sys.gradient);

%% Input checks

switch sys.gradient
    case 'xrmw', 
    case 'xrm',  
    case 'whole',
    case 'zoom',
    case 'hrmb',
    case 'hrmw',
    case 'magnus',
    otherwise, error('Gradient coil (%s) unkown', sys.gradient);
end


return

