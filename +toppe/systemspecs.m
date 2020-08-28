function sys = systemspecs(varargin)
% Create struct containing scanner hardware specs and related info.
%
% function sys = systemspecs(varargin)
%
% varargin is a comma-separated list of field/value pairs.
% See also struct2cellarg.m.
%
% Input options:
%   maxGrad         Max gradient amplitude. MUST match the physical gradient system limit on your scanner, 
%                   since gradient waveforms are scaled to this value. Default: 5 Gauss/cm
%   gradUnit        Gauss/cm (default) or mT/m 
%   maxSlew         Default: 20 Gauss/cm/ms. Can be less than the physical slew limit in your scanner.
%   slewUnit        Gauss/cm/ms (default) or T/m/s
%   maxRf           Default: 0.25 Gauss. Must be less than or equal to the maximum achievable RF
%                   amplitude for the transmit coil you're using.
%   rfUnit          Gauss (default) or mT 
%   raster          Default: 4e-6 sec
%   gamma           Default: 4.2576e3 Hz/Gauss
%   maxSlice        Max slice index in Pfile data storage. Don't yet know what the limit is here. Default: 200
%   maxView         Max view index in Pfile data storage. Also not sure here. Default: 500
%   maxEcho         Max echo index in Pfile data storage. Default: 16. 
%   addDelays       Set toppe.<start_core*/myrfdel/daqdel/timetrwait/timessi> = 0 (e.g., for converting to/from Pulseq)
%   toppe           struct describing the TOPPE interpreter/driver
%                      toppe.version          Default: 'v2'
%                      toppe.start_core_rf    Minimum start time (us) for rf modules. Default: 0. 
%                      toppe.start_core_daq   Minimum start time (us) for data acquisition modules. Default: 126
%                      toppe.start_core_grad  Minimum start time (us) for gradient-only modules. Default: 0.
%                      toppe.myrfdel          rf/gradient delay (us). Set to 'psd_rf_wait'.
%                      toppe.daqdel           daq/gradient delay (us). Set to 'psd_grd_wait'.
%                      toppe.timetrwait     Required delay at end of module (us). Determined empiricially. Default: 64.
%                      toppe.timessi        EPIC 'ssi' time, i.e., minimum duration/delay between modules (us). Default: 100.
%                      toppe.nMaxModules    max number of .mod files. Not known. Default: 30.
%                      toppe.nMaxWaveforms  max number of columns in waveform array. Not known. Default: 200.
%   gradient        Currently suppports 'xrmw', 'xrm', 'whole', 'zoom'. See toppe.pns(). Default: 'xrm'.
%
% Usage examples:
%  >> sys = systemspecs();                   % use default values
%  >> sys = system('maxSlice', 50);          % sets maxSlice to 50; otherwise contains default values

%% Defaults
maxGradDefault = 5;
sys.maxGrad  = maxGradDefault;
sys.gradUnit = 'Gauss/cm';
sys.maxSlew  = 20;
sys.slewUnit = 'Gauss/cm/ms';
maxRfDefault = 0.25;
sys.maxRf    = maxRfDefault;           % NB! Not sure what the hardware limit is here.
sys.rfUnit   = 'Gauss';
sys.raster   = 4e-6;           % sec
sys.gamma    = 4.2576e3;       % Hz/Gauss
sys.maxSlice = 2048;           % max dabslice. UI won't allow more than this to be entered
sys.maxView  = 600;            % also not sure here
sys.maxEcho  = 16;             % about right
sys.addDelays = true ;         % False: set time gaps to zero.

% sys.toppe struct relates to the TOPPE interpreter/driver.
% You probably shouldn't edit these.
sys.toppe.version = 'v2';
sys.toppe.start_core_rf  = 0;      % minimum start time (us) for rf modules
sys.toppe.start_core_daq = 126;    % minimum start time (us) for data acquisition modules
sys.toppe.start_core_grad = 0;     % minimum start time (us) for gradient-only modules
sys.toppe.myrfdel    = 78;         % rf/gradient delay (us) ( = psd_rf_wait). Inside scanner: 78. Outside scanner: 94.
sys.toppe.daqdel     = 84;         % daq/gradient delay (us) (= psd_grd_wait). Inside scanner: 84. Outside scanner: 100.
sys.toppe.timetrwait = 64;         % required delay at end of module (us)
sys.toppe.timessi    = 100;        % EPIC 'ssi' time, i.e., minimum duration/delay between modules
sys.toppe.nMaxModules   = 20;      % Not known at the moment, probably limited by total scanner system memory.
sys.toppe.nMaxWaveforms = 200;     % Not known at the moment.

% gradient subsystem. See also pns.m.
% Scanner  Gradient coil   chronaxie rheobase alpha  gmax  smax
% MR750w   XRMW            360d-6    20.0     0.324  33    120
% MR750    XRM             334d-6    23.4     0.333  50    200
% HDx      TRM WHOLE       370d-6    23.7     0.344  23    77
% HDx      TRM ZOOM        354d-6    29.1     0.309  40    150
%
% values on scanner from /w/config/Scandbdt.cfg + GRSubsystemHWO.xml
sys.gradient = 'xrm';

% Substitute specified system values as appropriate
sys = toppe.utils.vararg_pair(sys, varargin);

% Warn user if setting maxGrad or maxRf; this can lead to wrong gradient/b1 amplitudes on scanner.
if strcmp(sys.gradUnit, 'Gauss/cm')
	if sys.maxGrad ~= maxGradDefault
		warning(sprintf('Are you sure the gradients on your scanner are scaled to peak value %.1f G/cm?', sys.maxGrad));
	end
else
	if sys.maxGrad ~= 10*maxGradDefault
		warning(sprintf('Are you sure the gradients on your scanner are scaled to peak value %.1f mT/m?', sys.maxGrad));
	end
end
if sys.maxRf ~= maxRfDefault
	warning(sprintf('Are you sure you want to change maxRf? B1 scaling may be impacted. Recall that peak b1 is determined by maxRf value in tipdown.mod.'));
end

% If requested, set GE/EPIC-related time gaps to zero
if ~sys.addDelays
	sys.toppe.start_core = 0;
	sys.toppe.myrfdel    = 0;
	sys.toppe.daqdel     = 0;
	sys.toppe.start_core = 0;
	sys.toppe.timetrwait = 0;
	sys.toppe.timessi    = 0;
end

return

