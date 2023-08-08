function [rf, gx, gy, gz, tRange] = plotseq(sysGE, varargin)
% function [rf, gx, gy, gz, tRange] = plotseq(sysGE, varargin)
%
% Display pulse sequence, as specified in modules.txt and scanloop.txt in current folder
%
% Sequence timing here is approximate.
% See v6/figs/timing.svg in the TOPPEpsdSourceCode Github repo for detailed sequence timing.
%
% Inputs:
%   sysGE             struct specifying hardware system info, see systemspecs.m
%
% Options:
%   'timeRange'       Specify time range, [tStart tStop] (sec)
%   'blockRange'      Block range, [iStart iStop]
%   'doDisplay'       true/false
%
% Outputs:
%   tRange            [tStart tStop]  Start and stop times (sec)

%% parse inputs
arg.timeRange       = [];
arg.blockRange      = [];
arg.loop            = [];
arg.loopFile        = 'scanloop.txt';
arg.mods            = [];
arg.moduleListFile  = 'modules.txt';
arg.doDisplay       = true;
arg.doTimeOnly      = false;
arg.drawpause       = 1;
arg.gmax            = [];     % Gauss/cm
arg.rhomax          = [];  % Gauss
arg.printTime       = false;

arg = toppe.utils.vararg_pair(arg, varargin);

if isempty(arg.gmax)
    arg.gmax = sysGE.maxGrad;  % default
end
if isempty(arg.rhomax)
    arg.rhomax = sysGE.maxRF;  % default
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

%% Initialize counter and turn off display if we're only doing timings
if arg.doTimeOnly
    nsamples = 0;
    arg.doDisplay = false;
    for p = 1:size(modules,2) % Compute table of core durations
        core_size(p) = size(modules{p}.gx(:,1),1);
    end
end

%% Get block groups. Not used at the moment [TODO]
%blockGroups = toppe.readcoresfile('cores.txt');

%if arg.doTimeOnly 
%    gxlength = round(modules{p}.dur/raster);
%    nsamples = nsamples + gxlength;

%% build sequence. each sample is 4us.
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
            t = t + modules{p}.dur;
        end
    end
    tStart = t;
    for n = arg.blockRange(1):arg.blockRange(2)
        p = loop(n,1);   % module id
        if p == 0
            t = t + loop(n,14);    % us
        else
            t = t + modules{p}.dur;
        end
    end
    tStop = t;
    T = (tStart+raster/2) : raster : tStop;

    tRange = [tStart tStop]*1e-6;  % s

    [rho, th, gx, gy, gz, dur] = sub_getwavs(arg.blockRange(1), arg.blockRange(2), loop, modules, sysGE, arg.doTimeOnly);
else
    % Plot time range
    tic = arg.timeRange(1);
    toc = arg.timeRange(2);

    % find blocks containing start and stop times
    n = 1;
    t = 0;  % time counter (us)
    while t <= tic
        p = loop(n,1);   % module id
        i = loop(n,28);  % block group id

        if p == 0
            dur = loop(n,14);    % us
        else
            dur = modules{p}.dur;  % us
        end
        t = t + dur;
        n = n + 1;
    end

    tStart = max(t-dur, 0);
    nStart = max(n-1, 1);

    while t <= toc & n <= size(loop, 1)
        p = loop(n,1);   % module id
        i = loop(n,28);  % block group id

        if p == 0
            t = t + loop(n,14);    % us
        else
            t = t + modules{p}.dur;
        end
        n = n + 1;
    end

    tStop = t;
    nStop = max(n-1, 1);

    % get waveforms
    [rho, th, gx, gy, gz, dur] = sub_getwavs(nStart, nStop, loop, modules, sysGE, arg.doTimeOnly);

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
    end

    tRange = arg.timeRange*1e-6;
end

rf = rho.*exp(1i*th);


%if arg.printTime
%    fprintf(1, 'n %d: mindur = %d us, rf t = %d us, grad t = %d us\n', n, mindur, numel(rho)*raster, numel(gx)*raster);
%end

%% plot
if arg.doDisplay
%    T = (0:(numel(rho)-1))*raster/1000; % msec
    T = T*1e-3;  % ms
    
    Tend = max(T);
    
    gmax = 1.1*arg.gmax;   % Gauss/cm
    srho = 1.1*arg.rhomax; %max(1.1*max(abs(rho(:))),0.05);
    lw = 2;
    bgColor = 'k';

    t = tiledlayout(5, 1);
    ax1 = nexttile;
    plot(T, gx, '-y', 'LineWidth', lw);  ylabel('X (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax2 = nexttile;
    plot(T, gy, '-c', 'LineWidth', lw);  ylabel('Y (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);
    
    ax3 = nexttile;
    plot(T, gz, '-m', 'LineWidth', lw);  ylabel('Z (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax4 = nexttile;
    plot(T, rho, '-r', 'LineWidth', lw); ylabel('|b1| (Gauss)'); axis([T(1) Tend -srho srho]);
    set(gca, 'color', bgColor);  set(gca, 'XTick', []);

    ax5 = nexttile;
    plot(T, th, '-g', 'LineWidth', lw);  ylabel('∠b1 (rad)'); axis([T(1) Tend -1.1*pi 1.1*pi]);
    set(gca, 'color', bgColor);
    xlabel('time (ms)');

    t.TileSpacing = 'none';
    t.Padding = 'none';

    linkaxes([ax1 ax2 ax3 ax4 ax5], 'x');
end

% build sequence. Each sample is 4us.
function [rho, th, gx, gy, gz, dur] = sub_getwavs(blockStart, blockStop, loop, modules, sysGE, doTimeOnly)

rho = []; th = []; gx = []; gy = []; gz = [];
max_pg_iamp = 2^15-2;
raster = sysGE.raster;  % us

for n = blockStart : blockStop

    p = loop(n,1);   % module id
    i = loop(n,28);  % block group id
    w = loop(n,16);  % waveform index

    dur = 0;

    % Pure delay block
    if p == 0
        dur = dur + loop(n,14);  % TODO: make accurate
        if ~doTimeOnly   
            w = zeros(round(loop(n,14)/raster), 1);  % TODO: make accurate
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
    res = round(modules{p}.dur/raster);
    rho = [rho; rho1; zeros(res - modules{p}.res, 1)];
    th = [th; th1; zeros(res - modules{p}.res, 1)];
    gx = [gx; gx1; zeros(res - modules{p}.res, 1)];
    gy = [gy; gy1; zeros(res - modules{p}.res, 1)];
    gz = [gz; gz1; zeros(res - modules{p}.res, 1)];
end
