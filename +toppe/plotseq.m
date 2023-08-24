function [rf, gx, gy, gz, tRange] = plotseq(sysGE, varargin)
% function [rf, gx, gy, gz, tRange] = plotseq(sysGE, varargin)
%
% Display pulse sequence, as specified in modules.txt and scanloop.txt in current folder.
%
% Can also be called the old way for backward compatibility:
%     function plotseq(nStart, nStop, sysGE, varargin)
%
% NB! Sequence timing is approximate (for now).
% See v6/figs/timing.svg in the TOPPEpsdSourceCode Github repo for detailed sequence timing.
%
% Inputs:
%   sysGE             struct specifying hardware system info, see systemspecs.m
%
% Options:
%   'timeRange'       Specify time range, [tStart tStop] (sec)
%   'blockRange'      Block range, [iStart iStop]
%
% Outputs:
%   tRange            [tStart tStop]  Start and stop times (sec)
%
% Example:
%  Plot the whole sequence:
%    >> toppe.plotseq(sysGE, 'timeRange', [0 inf]);
%    or just:
%    >> toppe.plotseq(sysGE);

% In the code:
%   - time in us

%% parse input options

% Support the old way of calling plotseq: plotseq(nStart, nStep, sysGE, varargin)
if isnumeric(sysGE)
    arg.blockRange = [sysGE varargin{1}];
    sysGE = varargin{2};
    if length(varargin) > 2
        varargin = varargin(3:end);
    else
        varargin = {};
    end
else
    arg.blockRange = [];
end

% defaults
arg.timeRange       = [];
arg.loop            = [];
arg.loopFile        = 'scanloop.txt';
arg.mods            = [];
arg.moduleListFile  = 'modules.txt';
arg.doDisplay       = true;
arg.doTimeOnly      = false;
arg.gmax            = [];  % Gauss/cm
arg.rhomax          = [];  % Gauss
arg.printTime       = false;

arg = toppe.utils.vararg_pair(arg, varargin);

if isempty(arg.timeRange) & isempty(arg.blockRange)
    arg.timeRange = [0 inf];
end

arg.timeRange = round(arg.timeRange*1e6);   % us

if arg.doTimeOnly
    arg.doDisplay = false;
end

%% read scan files
% scanloop 
if isempty(arg.loop)
    loop = toppe.tryread(@toppe.readloop, arg.loopFile);
else
    loop = arg.loop;
end

% modules
if isempty(arg.mods)
    modules = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);
else
    modules = arg.mods;
end

% segment definitions (here referred to as 'module groups' in accordance with tv6.e)
if isfile('cores.txt')
    modGroups = toppe.tryread(@toppe.readcoresfile, 'cores.txt');
    useDefaultSegments = false;
else
    % Assign each .mod file to its own segment
    for p = 1:length(modules)
        modGroups{p} = p;
    end
    useDefaultSegments = true; 
end

%% Add end of segment label to loop so we know when to insert
%% the segmentRingdownTime (= 116 us)
n = 1;
isLastBlockInSegment = zeros(size(loop, 1), 1);  % initial values
while n < size(loop, 1)
    if useDefaultSegments
        i = loop(n, 1);   % segment ID = module ID
    else
        i = loop(n, end);   % segment ID
    end

    isLastBlockInSegment(n + length(modGroups{i}) - 1) = 1;
    n = n + length(modGroups{i});
end

%% Build sequence. Each time sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
raster = sysGE.raster;  % us

if ~isempty(arg.blockRange)
    % get start and stop times (for plotting)
    n = 1;
    t = 0;
    for n = 1:(arg.blockRange(1)-1)
        p = loop(n,1);   % module id
        if p == 0
            t = t + loop(n,14);    % us
        else
            minimumModuleDuration = modules{p}.res*raster;  % us
            t = t + max(modules{p}.dur, minimumModuleDuration);
        end
    end
    tStart = t;
    for n = arg.blockRange(1):arg.blockRange(2)
        p = loop(n,1);   % module id
        if p == 0
            t = t + loop(n,14);    % us
        else
            minimumModuleDuration = modules{p}.res*raster;  % us
            t = t + max(modules{p}.dur, minimumModuleDuration);
        end
    end
    tStop = t;
    T = (tStart+raster/2) : raster : tStop;

    tRange = [tStart tStop]*1e-6;  % s

    [rho, th, gx, gy, gz, dur] = sub_getwavs(arg.blockRange(1), arg.blockRange(2), loop, modules, sysGE, arg.doTimeOnly, isLastBlockInSegment, useDefaultSegments);
else
    % Plot time range
    tic = arg.timeRange(1);
    toc = arg.timeRange(2);

    % find blocks containing start and stop times
    n = 1;
    t = 0;  % time counter (us)
    dur = [];
    while t <= tic 
        if n == size(loop,1)-1
            error('invalid start time');
        end
        
        p = loop(n,1);   % module id
        if p == 0
            dur = loop(n,14);    % us
        else
            minimumModuleDuration = modules{p}.res*raster;  % us
            dur = max(modules{p}.dur, minimumModuleDuration);
        end
        t = t + dur;
        if isLastBlockInSegment(n)
            t = t + sysGE.segmentRingdownTime;
        end
        n = n + 1;
    end

    if isempty(dur)
        error('Invalid time range');
    end

    tStart = max(t-dur, 0);
    nStart = max(n-1, 1);

    while t < toc & n <= size(loop, 1)
        p = loop(n,1);   % module id

        if p == 0
            dur = loop(n,14);    % us
        else
            minimumModuleDuration = modules{p}.res*raster;  % us
            dur = max(modules{p}.dur, minimumModuleDuration);
        end
        t = t + dur;
        if isLastBlockInSegment(n)
            t = t + sysGE.segmentRingdownTime;
        end
        n = n + 1;
    end

    tStop = t;
    nStop = max(n-1, 1);

    % get waveforms
    [rho, th, gx, gy, gz, dur] = sub_getwavs(nStart, nStop, loop, modules, sysGE, arg.doTimeOnly, isLastBlockInSegment, useDefaultSegments);

    if ~arg.doTimeOnly
        % vector of time points (for plotting)
        % keep only desired time range
        T = (tStart+raster/2):raster:tStop;
        mask = T >= tic & T <= toc;
        T = T (mask);
        rho = rho(mask);
        th = th(mask);
        gx = gx(mask);
        gy = gy(mask);
        gz = gz(mask);
        tRange = [T(1) T(end)]*1e-6;
    else
        tRange = [tStart tStop]*1e-6;
    end

end

rf = rho.*exp(1i*th);

%% plot
if arg.doDisplay
    T = T*1e-3;  % ms
    
    % start/end of time axis (for plotting)
    tPad = (T(end)-T(1))/100;  % pad with some blank space before and after plot
    Tstart = T(1) - tPad;
    Tend = max(T) + tPad;

    if isempty(arg.gmax)
        gxmax = 1.3*max(abs(gx));
        if gxmax == 0
            gxmax = sysGE.maxGrad;
        end
        gymax = 1.3*max(abs(gy));
        if gymax == 0
            gymax = sysGE.maxGrad;
        end
        gzmax = 1.3*max(abs(gz));
        if gzmax == 0
            gzmax = sysGE.maxGrad;
        end
    else
        gxmax = arg.gmax;
        gymax = arg.gmax;
        gzmax = arg.gmax;
    end

    if isempty(arg.rhomax)
        rhomax = 1.3*max(abs(rho));
        if rhomax == 0
            rhomax = sysGE.maxRF;
        end
    else
        rhomax = arg.rhomax;
    end
    
    lw = 2;
    bgColor = 'k';

    t = tiledlayout(5, 1);
    ax1 = nexttile;
    plot(T, gx, '-y', 'LineWidth', lw);  ylabel('X (G/cm)'); axis([Tstart Tend -gxmax gxmax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax2 = nexttile;
    plot(T, gy, '-c', 'LineWidth', lw);  ylabel('Y (G/cm)'); axis([Tstart Tend -gymax gymax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);
    
    ax3 = nexttile;
    plot(T, gz, '-m', 'LineWidth', lw);  ylabel('Z (G/cm)'); axis([Tstart Tend -gzmax gzmax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax4 = nexttile;
    plot(T, rho, '-r', 'LineWidth', lw); ylabel('|b1| (Gauss)'); axis([Tstart Tend -rhomax rhomax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax5 = nexttile;
    plot(T, th, '-g', 'LineWidth', lw);  ylabel('∠b1 (rad)'); axis([T(1) Tend -1.1*pi 1.1*pi]);
    set(gca, 'color', bgColor);
    xlabel('time (ms)');

    t.TileSpacing = 'none';
    t.Padding = 'none';

    linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');  % common zoom setting (along time axis) for all tiles
end


%% Function to build sequence. Each sample is 4us.
function [rho, th, gx, gy, gz, dur] = sub_getwavs(blockStart, blockStop, loop, modules, sysGE, doTimeOnly, isLastBlockInSegment, useDefaultSegments)

rho = []; th = []; gx = []; gy = []; gz = [];
max_pg_iamp = 2^15-2;
raster = sysGE.raster;  % us

for n = blockStart : blockStop

    p = loop(n,1);   % module id
    w = loop(n,16);  % waveform index

    dur = 0;

    % Pure delay block
    if p == 0
        dur = dur + loop(n,14);
        if isLastBlockInSegment(n)
            dur = dur + sysGE.segmentRingdownTime;
        end
        if ~doTimeOnly   
            w = zeros(round(dur/raster), 1);
            rho = [rho; w];
            th = [th; w];
            gx = [gx; w];
            gy = [gy; w];
            gz = [gz; w];
        end
        continue;
    end

    if modules{p}.hasRF
        ia_rf = loop(n,2);
    else
        ia_rf = 0;
    end
    ia_th = loop(n,3);
    ia_gx = loop(n,4);
    ia_gy = loop(n,5);
    ia_gz = loop(n,6);

    % Calculate RF and gradient waveforms for this block
    rho1 = [ia_rf/max_pg_iamp*abs(modules{p}.rf(:,w))];
    th1 = [ia_rf/max_pg_iamp*angle(modules{p}.rf(:,w))];
    gxit = ia_gx/max_pg_iamp*modules{p}.gx(:,w);
    gyit = ia_gy/max_pg_iamp*modules{p}.gy(:,w);
    gzit = ia_gz/max_pg_iamp*modules{p}.gz(:,w);

    % apply 3d rotation matrix 
    % (which also accounts for any in-plane 2D rotation, i.e., 'phi' in write2loop.m)
    Rv = loop(n,17:25)/max_pg_iamp;  % stored in row-major order
    R = reshape(Rv, 3, 3);
    G = R * [gxit(:)'; gyit(:)'; gzit(:)'];
    gx1 = G(1,:)';
    gy1 = G(2,:)';
    gz1 = G(3,:)';

    % apply RF phase offset
    if modules{p}.hasRF
        th1 = th1 + loop(n,12)/max_pg_iamp*pi;
        th1 = angle(exp(1i*th1));   % wrap to [-pi pi] range
    end

    % pad to desired duration and add to running waveform
    minimumModuleDuration = modules{p}.res*raster;  % us
    dur = max(modules{p}.dur, minimumModuleDuration);
    if isLastBlockInSegment(n)
        dur = dur + sysGE.segmentRingdownTime;
    end
    res = round(dur/raster);
    rho = [rho; rho1; zeros(res - modules{p}.res, 1)];
    th = [th; th1; zeros(res - modules{p}.res, 1)];
    gx = [gx; gx1; zeros(res - modules{p}.res, 1)];
    gy = [gy; gy1; zeros(res - modules{p}.res, 1)];
    gz = [gz; gz1; zeros(res - modules{p}.res, 1)];
end
