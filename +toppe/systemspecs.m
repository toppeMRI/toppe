function sys = systemspecs(varargin)
% Create struct containing scanner hardware specs and related info.
%
% function sys = systemspecs(varargin)
% 
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
%
% Times are in microseconds


%% Defaults

% Global constants
sys.raster = 4;                 % us. Raster time for gradient and RF waveforms.
sys.gamma  = 42.576e6;          % Hz/T

% Scanner-specific settings
sys.B0 = 3.0;                     % field strength (T)
sys.gradient = 'xrm';             % gradient coil
sys.psd_rf_wait = 148;            % rf/gradient delay (us)
sys.psd_grd_wait = 156;           % ADC/gradient delay (us).
%sys.segmentDeadTime = 12;         % Dead time before start of block group, equals RUP_GRD(9us)
sys.segmentRingdownTime = 116;    % Delay at end of block group, equals 4us + timssi. 
sys.forbiddenEspRange = [410 510];    % (us) Forbidden echo spacings (mechanical resonance). See /usr/g/bin/epiesp.dat
sys.tminwait   = 12;              % minimum duration of wait pulse in EPIC code (us)

% Design choices (need not equal scanner limits)
sys.maxGrad  = 4;      % Gauss/cm
sys.maxSlew  = 15;     % Gauss/cm/ms
sys.maxRF    = 0.25; 
sys.rfDeadTime = 72;         % us. Must be >= 72us
sys.rfRingdownTime = 54;     % us. Must be >= 54us
sys.adcDeadTime = 40;        % us. Must be >= 40us

% The following determine the slice/echo/view indexing in the data file
sys.maxSlice = 2048;           % max dabslice. UI won't allow more than this to be entered
sys.maxView  = 2048;
sys.maxEcho  = 1;             % actual value seems to be 16, but we won't use this dimension so only allow 1 here


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

