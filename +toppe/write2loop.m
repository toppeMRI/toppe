function write2loop(modname, system, varargin)
% function write2loop(modname, system, varargin)
%
% Writes a module and settings to the scanloop file for a sequence.
% To create a scanloop, first initialize the file using:
%    write2loop('setup')
%
% Then write RF or DAQ modules in order they will be played:
%   Ex.   write2loop('tipdown.mod',sys,'RFspoil',true,'RFoffset',-5);
%         write2loop('readout.mod',sys,'slice',1,'view',80,'echo',5);
%
% Complete the scanloop using:
%    write2loop('finish')
%
%    Without calling 'finish', your sequence will not run!
%
% Inputs:
%    modname            [string] Name of module to call, ex. "tipdown.mod"
%                                OR initialization/completion calls
%                                "setup" / "finish"
%    system             [1 1]    struct specifying hardware system info, see systemspecs.m
% Regular input options:
%    loopFile           [string] Loop file to write to, default: 'scanloop.txt'
%    moduleListFile     [string] File that lists modules to use, default: 'modules.txt'
%    Gamplitude         [3 1]    Amplitude scaling of gradient waveform
%                                Range -1 to 1, default: 1
%    core               [1 1]    core id (positive int; defined in cores.txt)
%                                Default: 1
%    waveform           [1 1]    Waveform number to play (1,...,nwaveforms)
%                                Default: 1
%    textra             [1 1]    Extra time at end of waveform (ms)
%    rot                [1 1]    In-plane gradient rotation angle (radians), around the axis of rotation defined by rotmat.
%    rotmat             [3 3]    3x3 rotation matrix (for toppe >= v3). If 'version'=2, then 'rot' is applied and 'rotmat' is ignored.
%    trig               [1 1]    Trigger mode (int): internal (0) or cardiac (1)
%    trigout            [1 1]    Trigger out (int): off (0) or on (1)
%    version            [1 1]    Default: 4. Only needs to be passed if the first argument is 'setup' (otherwise ignored w/o warning)
%
% RF module input options:
%    RFamplitude        [1 1]    Amplitude scaling of RF waveform
%                                Range -1 to 1, default: 1
%    RFphase            [1 1]    Phase of RF (rad), default 0
%    RFspoil            [1 1]    Automatically use and update RF spoiling 
%                                phase, using an angle of 117 deg.
%                                True/false, default false
%    RFoffset           [1 1]    Integer frequency offset for RF pulse (Hz)
%
% DAQ module input options:
%    slice              [1 1]    Slice number to write data to
%                                Range 1-2048, or 'dis' to discard (default 1)
%    echo               [1 1]    Echo number to write data to
%                                Range 1-??? (default 1)
%    view               [1 1]    View number to write data to
%                                Range 1-??? (default 1)
%    dabmode            [string] Acquisition mode. Options:
%                                'on'    Regular data acq (default)
%                                'off'   No acq
%                                'add'   Average with existing
%                                'reset' Erase data (?)

import toppe.*
import toppe.utils.*

%% parse inputs
% Defaults
arg.loopFile        = 'scanloop.txt';
arg.moduleListFile  = 'modules.txt';
arg.Gamplitude      = [1 1 1]';
arg.core            = 1;
arg.waveform        = 1;
arg.textra          = 0;
arg.RFamplitude     = 1;
arg.RFphase         = 0;
arg.DAQphase        = 0;
arg.RFspoil         = false;
arg.RFoffset        = 0;
arg.slice           = 1;
arg.echo            = 1;
arg.view            = 1;
arg.dabmode         = 'on';
arg.rot             = 0;
arg.rotmat          = eye(3);
arg.trig            = trig_intern;
arg.version         = 4;

% Adding trigger out for toppe v5
arg.trigout         = 0;

% substitute with explicit inputs
arg = toppe.utils.vararg_pair(arg, varargin);

% check inputs
checkInputs(arg);

%% Apply in-plane rotation angle 'rot' to arg.rotmat
n = arg.rotmat*[0 0 1]';   % axis of rotation
R = angleaxis2rotmat(arg.rot, n);
rotmat = R*arg.rotmat;

%% Set up persistant values for tracking values inbetween calls
persistent rf_spoil_seed_cnt
persistent rf_spoil_phase
persistent irfphase
persistent idaqphase
persistent setupdone
persistent modules
persistent d
persistent d_index
persistent toppeVer

%% If modname is setup, then init the file
if strcmp(modname,'setup')    

    toppeVer = arg.version;

    % Zero out spoiling angle tracking variables
    rf_spoil_phase = 0;
    irfphase = 0;
    rf_spoil_seed_cnt = 0;
    
    % Load modules in only once and keep persistent for speed
    modules = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);

    % Load cores
    
    % Preallocate 
    switch toppeVer
        case 6
            d = zeros(1000000,28); % 28th column is core id
        case 5
            d = zeros(1000000,27); % 27th column is trigger out flag (0 or 1)
        case 4
            d = zeros(1000000,26);   % 26th column is trigger in mode (0 or 1)
        case 3
            d = zeros(1000000,25);   % Columns 17-25 contain 3D rotation matrix
        otherwise
            d = zeros(1000000,16);
    end
    d_index = 1; % Current line to write to
    
    % Mark we're set up
    setupdone = 1;
    return
end

%if strcmp(modname, 'core')

if strcmp(modname,'finish')
    % Remove rows that are all zeros (preallocated but unused)
    d = d(1:d_index-1,:);
    
    % Check if all lines are integers
    if ~all(d == round(d))
        error('Value in d is a non-integer, this should never happen!')
    end
    
    % Calculate header values
    nt = size(d,1);              % number of startseq() calls
    maxslice = max(d(:,7));
    maxecho = max(d(:,8));
    maxview = max(d(:,9));

    % check if max 'slice', 'echo', and 'view' numbers in scanloop.txt exceed system limits
    if maxslice > system.maxSlice
        warning('maxslice > system.maxSlice -- scan may not run!'); 
    end
    if maxecho + 1 > system.maxEcho    % +1 since 'echo' starts at 0
        warning('maxecho > system.maxEcho -- scan may not run!');   
    end
    if maxview > system.maxView
        warning('maxview > system.maxView -- scan may not run!');   
    end

    dur = toppe.getscantime(system,'loopArr',d,'mods',modules);
    udur = round(dur * 1e6);
    if toppeVer > 2
        newParams = [nt maxslice maxecho maxview udur toppeVer];
    else
        newParams = [nt maxslice maxecho maxview udur];
    end
    
    % Write file header
    fid = fopen(arg.loopFile, 'w+', 'ieee-be');
    if toppeVer > 2
        fprintf(fid, 'nt\tmaxslice\tmaxecho\tmaxview\tscandur\tversion\n');
    else
        fprintf(fid, 'nt\tmaxslice\tmaxecho\tmaxview\tscandur\n');
    end
    if toppeVer > 2
        lineformat = '%d\t%d\t%d\t%d\t%d\t%d\n';
    else
        lineformat = '%d\t%d\t%d\t%d\t%d\n';
    end
    fprintf(fid, lineformat, newParams); % Updated line
    fprintf(fid, headerline);
    fclose(fid);

    % Write scanloop lines from d matrix.
    % Precision=8 needed to avoid large numbers written in scientific notation.
    dlmwrite(arg.loopFile, d, '-append','delimiter', '\t', 'precision', 8);  
    
    % Disable setup flag to conclude file
    modules = []; % Erase modules struct from mem for safety in case .mod files get updated
    setupdone = [];
    return
end

%% Check for setup 
if isempty(setupdone)
    error('Set up not performed, call write2loop(''setup'') before writing module files.');
end

%% Find input module in persistant module struct
% Find module in cell array and store the index in moduleno
if strcmp(modname, 'delay')
    % pure delay block
    iModule = 0;
    module = struct();
    module.hasRF = 0;
    module.hasDAQ = 0;
else
    for iModule = 1:size(modules,2)
        if strcmp(modname,modules{iModule}.fname)
            break;
        end
        if iModule == size(modules,2)
            error(['Can''t find ' modname ' in ' arg.moduleListFile]);
        end
    end
    module = modules{iModule};
end

% Calculate gradient scalings
ia_gx = 2*round(arg.Gamplitude(1)*max_pg_iamp/2);
ia_gy = 2*round(arg.Gamplitude(2)*max_pg_iamp/2);
ia_gz = 2*round(arg.Gamplitude(3)*max_pg_iamp/2);

% Calculate time at end of module
textra_us = 4*round(1000*arg.textra/4);

% 3d rotation matrix
if toppeVer > 2
    if ndims(rotmat) ~= 2 | ~all(size(rotmat)==3) % | norm(rotmat) ~= 1
        error('rotmat must be a 3x3 orthonormal matrix');
    end
    rt = rotmat';
    drot = round(max_pg_iamp*rt(:)');   % vectorized, in row-major order
    phi = 0; % in-plane rotation (already applied to 3D rotation matrix)
else
    drot = [];
    phi = angle(exp(1i*arg.rot));     % inplane rotation. wrap to [-pi pi]
end;

% 2D in-plane rotation ('iphi' only applies to v2)
iphi = 2*round(phi/pi*max_pg_iamp/2);

% trigger
if toppeVer > 3
    trig = arg.trig;
else
    trig = []; 
end

% trigger (TTL pulse) out
if toppeVer > 4
    trigout = arg.trigout;
end

if module.hasRF % Write RF module
    % Do RF amplitude stuff
    ia_rf = 2*round(arg.RFamplitude*max_pg_iamp/2);
    ia_th = max_pg_iamp;
    
    % Do RF phase stuff
    if arg.RFspoil % If spoil is called, replace RF phase with spoil phase
        updateSpoilPhase();
        irfphase = phase2int(rf_spoil_phase);
        idaqphase = irfphase;
    else
        irfphase = phase2int(arg.RFphase);
        idaqphase = phase2int(arg.DAQphase);
    end
    
    % Dummy slice/echo/view
    dabslice = 0;
    dabecho = 0;
    dabview = 1;
    
    % Frequency offset in Hz
    if arg.RFoffset ~= round(arg.RFoffset)
        warning('Rounding frequency offset to nearest value');
        arg.RFoffset = round(arg.RFoffset);
    end
    f = arg.RFoffset;
 
    % Write line values and increment
    switch toppeVer
        case {4}
            d(d_index,:) = [iModule ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview 0 iphi irfphase irfphase textra_us f arg.waveform drot trig];
        case {5}   % Adding 1 column for trigger out.
            d(d_index,:) = [iModule ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview 0 iphi irfphase irfphase textra_us f arg.waveform drot trig trigout];
        case {6}   % Adding 1 column for core id
            d(d_index,:) = [iModule ia_rf ia_th ia_gx ia_gy ia_gz dabslice dabecho dabview 0 iphi irfphase irfphase textra_us f arg.waveform drot trig trigout arg.core];
    end    
    d_index = d_index + 1;
elseif module.hasDAQ % Write DAQ module
    % Set slice/echo/view
    if all(arg.slice == 'dis')
        dabslice = 0;
    elseif arg.slice > 0
        dabslice = arg.slice;
    else
        error('Slice must be ''dis'' or > 0');
    end
    dabecho = arg.echo-1; % Index echo from 0 to n-1
    dabview = arg.view;
    % Xiaoxi: update for adc freq offset
    f = arg.RFoffset;
    
    % receive phase
    if arg.RFspoil % If spoil is called, replace RF phase with spoil phase
       idaqphase = irfphase;
    else
       idaqphase = phase2int(arg.DAQphase);
    end
    % Xiaoxi: save adc freq offset in rf freq offset column
    switch toppeVer
        case {4} 
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz dabslice dabecho dabview dabval(arg.dabmode) iphi idaqphase idaqphase textra_us f arg.waveform drot trig];
        case {5} 
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz dabslice dabecho dabview dabval(arg.dabmode) iphi idaqphase idaqphase textra_us f arg.waveform drot trig trigout];
        case {6} 
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz dabslice dabecho dabview dabval(arg.dabmode) iphi idaqphase idaqphase textra_us f arg.waveform drot trig trigout arg.core];
    end    

    d_index = d_index + 1;
else
    % gradients only
    %error(['Module didn''t have RF or DAQ specified, check ' arg.moduleListFile]);

    % rotation
    switch toppeVer
        case {4}
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz 0 0 0 0 iphi 0 0 textra_us 0 arg.waveform drot trig];
        case {5}
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz 0 0 0 0 iphi 0 0 textra_us 0 arg.waveform drot trig trigout];
        case {6}
            d(d_index,:) = [iModule 0 0 ia_gx ia_gy ia_gz 0 0 0 0 iphi 0 0 textra_us 0 arg.waveform drot trig trigout arg.core];
    end
    d_index = d_index + 1;
end
return

%% Nested functions
    function iphs = phase2int(phs)
        phstmp = atan2(sin(phs), cos(phs));          % wrap phase to (-pi,pi) range
        iphs = 2*round(phstmp/pi*max_pg_iamp/2);     % short int
    end

    function updateSpoilPhase
        rf_spoil_seed_cnt = rf_spoil_seed_cnt + 1;
        rf_spoil_phase = rf_spoil_phase + (rf_spoil_seed/180 * pi)*rf_spoil_seed_cnt;  % radians
    end

end

%% Other functions
function checkInputs(arg)
% Check for valid inputs to function arguments
if ~ischar(arg.loopFile) || isempty(arg.loopFile)
    error('Loop file argument must be a non-empty character array');
end
if ~ischar(arg.moduleListFile) || isempty(arg.moduleListFile)
    error('Module list argument must be a non-empty character array');
end
if any(abs(arg.Gamplitude)>1)
    error('Gradient amplitudes must be specified as a value between -1 and 1');
end
if any(size(arg.Gamplitude) ~= [3 1])
    error('Gradient amplitudes must be specified as a 3x1 vector');
end
if arg.waveform ~= round(abs(arg.waveform))
    error('Waveform must be specified as an integer > 1');
end
if arg.textra < 0
    error('textra must be non-negative');
end
if abs(arg.RFamplitude) > 1
    error('RF amplitude must be specified as a value between -1 and 1');
end
if ~(arg.RFspoil == 0 || arg.RFspoil == 1)
    error('RF spoil flag must be T/F');
end
if arg.RFspoil == 1 && arg.RFphase ~= 0
    error('RFspoil and RFphase were both specified, only use one at a time');
end
if arg.RFoffset ~= round(arg.RFoffset)
    error('RF Offset must be specified as an integer');
end
if ischar(arg.slice)
    if ~strcmp(arg.slice,'dis')
        error('Slice must be specified as an integer > 1 or ''dis''');
    end
elseif any(arg.slice ~= round(arg.slice))
    error('Slice must be an integer > 0');
end
if arg.echo ~= round(abs(arg.echo)) || arg.echo < 1
    error('Echo must be an integer > 0');
end
if arg.view ~= round(abs(arg.view)) || arg.view < 1
    error('View must be an integer > 0');
end
if arg.version < 4 || arg.version > 6
    error('version must be 4, 5, or 6');
end
end

%% Define constants via functions
function h = max_pg_iamp
h = 2^15-2; %Maximum single int gradient value
end

function h = rf_spoil_seed
h = 117;
end

function t = trig_intern
t = 0;
end

function t = trig_ecg
t = 1;
end

function h = dformat
%h = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n';
h = '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d \n';
end

function h = headerline
h = 'Core a_rf a_th a_gx a_gy a_gz dabslice dabecho dabview dabmode rot rfphase recphase textra freq waveform rotmat...\n';
end

function h = dabval(dabmode)
switch dabmode
    case 'off'
        h = 0;
    case 'on'
        h = 1;
    case 'add'
        h = 2;
    case 'reset'
        h = 3;
    otherwise
        error('Invalid dabmode option specified');
end
end
