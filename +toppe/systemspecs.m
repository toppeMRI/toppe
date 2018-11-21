function sys = systemspecs(varargin)
% Create struct containing scanner hardware specs and related info.
%
% function sys = systemspecs(varargin)
%
% varargin is a comma-separated list of field/value pairs.
% See also struct2cellarg.m.
%
% Input options:
%   maxGrad         Default: 5 Gauss/cm
%   gradUnit        mT/m (default) or Gauss/cm. 
%   maxSlew         20 Gauss/cm/ms
%   slewUnit        Gauss/cm/ms (default) or T/m/s
%   maxRf           Default: 0.25 Gauss
%   rfUnit          Gauss (default) or mT 
%   raster          Default: 4e-6 sec
%   gamma           Default: 4.2576e3 Hz/Gauss
%   maxSlice        Max slice index in Pfile data storage. I don't yet know what the limit is here. Default: 200
%   maxView         Max view index in Pfile data storage. Also not sure here. Default: 500
%   maxEcho         Max echo index in Pfile data storage. Default: 16. 
%   addDelays       Set toppe.<start_core/myrfdel/daqdel/timetrwai/timessi> = 0 (e.g., for converting to/from Pulseq)
%   toppe           struct describing the TOPPE interpreter/driver
%                      toppe.version        Default: 'v2'
%                      toppe.start_core     Minimum rf/gradient start time (us). Default: 224.
%                      toppe.myrfdel        rf/gradient delay (us). Default: 94.
%                      toppe.daqdel         daq/gradient delay (us). Default: 80.
%                      toppe.timetrwait     Required delay at end of module (us). Default: 64.
%                      toppe.timessi        EPIC 'ssi' time, i.e., minimum duration/delay between modules (us). Default: 100.
%                      toppe.nMaxModules    max number of .mod files. Not known. Default: 30.
%                      toppe.nMaxWaveforms  max number of columns in waveform array. Not known. Default: 200.
%
% Usage examples:
%  >> sys = systemspecs();                   % use default values
%  >> sys = system('maxSlice', 50);          % sets maxSlice to 50; otherwise contains default values
%
% $Id: systemspecs.m,v 1.10 2018/11/12 20:20:33 jfnielse Exp $
% $Source: /export/home/jfnielse/Private/cvs/projects/psd/toppe/matlab/+toppe/systemspecs.m,v $

%% Defaults
sys.maxGrad  = 5;
sys.gradUnit = 'Gauss/cm';
sys.maxSlew  = 20;
sys.slewUnit = 'Gauss/cm/ms';
sys.maxRf    = 0.25;           % NB! Not sure what the hardware limit is here.
sys.rfUnit   = 'Gauss';
sys.raster   = 4e-6;           % sec
sys.gamma    = 4.2576e3;       % Hz/Gauss
sys.maxSlice = 200;            % max dabslice. I don't yet know what the limit is here.
sys.maxView  = 250;            % also not sure here
sys.maxEcho  = 16;             % about right
sys.addDelays = true ;         % False: set time gaps to zero.

% sys.toppe struct relates to the TOPPE interpreter/driver.
% You shouldn't edit these.
sys.toppe.version = 'v2';
sys.toppe.start_core = 224;        % minimum rf/gradient start time (us)
sys.toppe.myrfdel    = 94;         % rf/gradient delay (us) ( = psd_rf_wait)
sys.toppe.daqdel     = 100;        % daq/gradient delay (us) (= psd_grd_wait)
sys.toppe.timetrwait = 64;         % required delay at end of module (us)
sys.toppe.timessi    = 100;        % EPIC 'ssi' time, i.e., minimum duration/delay between modules
sys.toppe.nMaxModules   = 20;      % Not known at the moment, probably limited by total scanner system memory.
sys.toppe.nMaxWaveforms = 200;     % Not known at the moment.

% Substitute specified system values as appropriate
sys = toppe.utils.vararg_pair(sys, varargin);

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

