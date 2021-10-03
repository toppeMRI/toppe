function sys = systemspecs(varargin)
% Create struct containing scanner hardware specs and related info.
%
% function sys = systemspecs(varargin)
% 
% 'maxSlew' and 'maxGrad' options can be < scanner limit, and can vary across .mod files. 
% The usual mode of operation in TOPPE is to use the default value for 
% maxRF (0.25G) for all .mod files, to ensure accurate b1 scaling in the interpreter.
%
% Input options:
%   maxGrad         Max gradient amplitude. Can be less than the physical limit for your scanner. 
%   gradUnit        Gauss/cm (default) or mT/m 
%   maxSlew         Default: 20 Gauss/cm/ms. Can be less than the physical slew limit for your scanner.
%   slewUnit        Gauss/cm/ms (default) or T/m/s
%   maxRF           Default: 0.25 Gauss. Using the default should ensure correct b1 scaling.
%   rfUnit          Gauss (default) or mT 
%   raster          Default: 4e-6 sec
%   gamma           Default: 4.2576e3 Hz/Gauss
%   maxSlice        Max slice index in Pfile data storage. Don't yet know what the limit is here. Default: 200
%   maxView         Max view index in Pfile data storage. Also not sure here. Default: 500
%   maxEcho         Max echo index in Pfile data storage. Default: 16. 
%   addDelays       Set toppe.<start_core*/myrfdel/daqdel/timetrwait/timessi> = 0 (e.g., for converting to/from Pulseq)
%   version         Default: 'v4'
%   gradient        For PNS calculation (see toppe.pns()). Currently suppports:
%                   'xrm' (MR750) (default)
%                   'xrmw' (MR750W) 
%                   'whole' (HDx) 
%                   'zoom' (HDx) 
%                   'hrmb' (UHP)
%   start_core_rf   Minimum start time (us) for rf modules. Default: 0. 
%   start_core_daq  Minimum start time (us) for data acquisition modules. Default: 126
%   start_core_grad Minimum start time (us) for gradient-only modules. Default: 0.
%   myrfdel         rf/gradient delay (us). Set to 'psd_rf_wait'.
%   daqdel          adc/gradient delay (us). Set to 'psd_grd_wait'.
%   timetrwait      Required delay at end of module (us). Determined empiricially. Default: 64.
%   timessi         EPIC 'ssi' time, i.e., minimum duration/delay between modules (us). Default: 100.
%   nMaxModules     max number of .mod files. Not known. Default: 30.
%   nMaxWaveforms   max number of columns in waveform array. Not known. Default: 200.
%
% Example 1 (typical usage):
% >> GE.sys = toppe.systemspecs('maxSlew', 12.3, 'slewUnit', 'Gauss/cm/ms', ...
% >> 'maxGrad', 5, 'gradUnit', 'Gauss/cm');
% 
% Example 2: Specify detailed timing info to match scanner:
% >> GE.sys = toppe.systemspecs('maxSlew', 12.3, ...
%    'maxGrad', 5, ...
%    'start_core_daq', 126, ...
%    'myrfdel', 94, ...
%    'daqdel', 100, ...
%    'timetrwait', 0, ...
%    'timessi', 100);


%% Defaults

% Hardware/design specs
maxRFDefault = 0.25;
sys.maxGrad  = 4;
sys.gradUnit = 'Gauss/cm';
sys.maxSlew  = 10;
sys.slewUnit = 'Gauss/cm/ms';
sys.maxRF    = maxRFDefault;           % NB! Not sure what the hardware limit is here.
sys.rfUnit   = 'Gauss';
sys.raster   = 4e-6;           % sec
sys.gamma    = 4.2576e3;       % Hz/Gauss
sys.maxSlice = 2048;           % max dabslice. UI won't allow more than this to be entered
sys.maxView  = 600;            % also not sure here
sys.maxEcho  = 16;             % about right
sys.addDelays = true ;         % False: set time gaps to zero.

% TOPPE version
sys.version = 'v4';

% timing
sys.start_core_rf  = 0;      % minimum start time (us) for rf modules
sys.start_core_daq = 126;    % minimum start time (us) for data acquisition modules
sys.start_core_grad = 0;     % minimum start time (us) for gradient-only modules
sys.myrfdel    = 78;         % rf/gradient delay (us) ( = psd_rf_wait). Inside scanner: 78. Outside scanner: 94.
sys.daqdel     = 84;         % daq/gradient delay (us) (= psd_grd_wait). Inside scanner: 84. Outside scanner: 100.
sys.timetrwait = 64;         % required delay at end of module (us)
sys.timessi    = 100;        % EPIC 'ssi' time, i.e., minimum duration/delay between modules
sys.tminwait   = 12;         % minimum duration of wait pulse in EPIC code

% gradient subsystem. See also pns.m.
% Scanner  Gradient coil   chronaxie rheobase alpha  gmax  smax
% MR750w   XRMW            360d-6    20.0     0.324  33    120
% MR750    XRM             334d-6    23.4     0.333  50    200
% HDx      TRM WHOLE       370d-6    23.7     0.344  23    77
% HDx      TRM ZOOM        354d-6    29.1     0.309  40    150
%
% values on scanner from /w/config/Scandbdt.cfg + GRSubsystemHWO.xml
sys.gradient = 'xrm';

% other
sys.nMaxModules   = 20;      % Not known at the moment, probably limited by total scanner system memory.
sys.nMaxWaveforms = 200;     % Not known at the moment.

%% Substitute specified system values
sys = toppe.utils.vararg_pair(sys, varargin);

if sys.maxRF ~= maxRFDefault
    warning('Using non-default max RF setting -- b1 scaling may be incorrect!');
end

% If requested, set GE/EPIC-related time gaps to zero
if ~sys.addDelays
	sys.start_core_rf = 0;
	sys.start_core_daq = 0;
	sys.start_core_grad = 0;
	sys.myrfdel    = 0;
	sys.daqdel     = 0;
	sys.timetrwait = 0;
	sys.timessi    = 0;
end

return

