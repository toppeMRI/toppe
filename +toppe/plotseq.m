function [rf, gx, gy, gz] = plotseq(nstart, nstop, sysGE, varargin)
% function [rf, gx, gy, gz] = plotseq(nstart, nstop, sysGE, varargin)
%
% Display pulse sequence, as specified in modules.txt, scanloop.txt, and timing.txt.
%
% Sequence timing here is approximate.
% See v6/figs/timing.svg in the TOPPEpsdSourceCode Github repo for detailed sequence timing.
%
% Inputs:
%   nstart,nstop       first and last startseq calls (as specified in scanloop.txt)
%   sysGE             struct specifying hardware system info, see systemspecs.m

%% parse inputs
arg.loopArr         = [];
arg.loopFile        = 'scanloop.txt';
arg.mods            = [];
arg.moduleListFile  = 'modules.txt';
arg.doDisplay       = true;
arg.doTimeOnly      = false;
arg.drawpause       = 1;
arg.gmax            = 5;     % Gauss/cm
arg.rhomax          = 0.25;  % Gauss
arg.printTime       = false;

arg = toppe.utils.vararg_pair(arg, varargin);

%% read scan files as needed
% scanloop array
if isempty(arg.loopArr)
    loopArr = toppe.tryread(@toppe.readloop, arg.loopFile);
else
    loopArr = arg.loopArr;
end

% module waveforms
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

%% Get block groups
blockGroups = toppe.readcoresfile('cores.txt');

%% build sequence. each sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
dt = 4;  % us
max_pg_iamp = 2^15-2;
raster = sysGE.raster;  % us

for n = nstart:nstop

    i = loopArr(n,28);  % block group id
    p = loopArr(n,1);   % module id
    w = loopArr(n,16);  % waveform index

    % Pure delay block
    if p == 0
       w = zeros(round(loopArr(n,14)/raster), 1);  % TODO: make accurate
       rho = [rho; w];
       th = [th; w];
       gx = [gx; w];
       gy = [gy; w];
       gz = [gz; w];
       continue;
    end

    if modules{p}.hasRF
        ia_rf = loopArr(n,2);
    else
        ia_rf = 0;
    end
    ia_th = loopArr(n,3);
    ia_gx = loopArr(n,4);
    ia_gy = loopArr(n,5);
    ia_gz = loopArr(n,6);

    if arg.doTimeOnly % Calculate the length of one waveform and add it to our sample counter
        gxlength = round(modules{p}.dur/raster);
        nsamples = nsamples + gxlength;
    else % Calculate RF and gradients as normal
        % waveforms for this block
        rho1 = [ia_rf/max_pg_iamp*abs(modules{p}.rf(:,w))];
        th1 = [ia_rf/max_pg_iamp*angle(modules{p}.rf(:,w))];
        gxit = ia_gx/max_pg_iamp*modules{p}.gx(:,w);
        gyit = ia_gy/max_pg_iamp*modules{p}.gy(:,w);
        gzit = ia_gz/max_pg_iamp*modules{p}.gz(:,w);

        % apply 3d rotation matrix 
        % (which also accounts for any in-plane 2D rotation, i.e., 'phi' in write2loop.m)
        Rv = loopArr(n,17:25)/max_pg_iamp;  % stored in row-major order
        R = reshape(Rv, 3, 3);
        G = R * [gxit(:)'; gyit(:)'; gzit(:)'];
        gx1 = G(1,:)';
        gy1 = G(2,:)';
        gz1 = G(3,:)';

        % apply RF phase offset
        if modules{p}.hasRF
            th1 = th1 + loopArr(n,12)/max_pg_iamp*pi;
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
    
    if arg.printTime
        fprintf(1, 'n %d: mindur = %d us, rf t = %d us, grad t = %d us\n', n, mindur, numel(rho)*raster, numel(gx)*raster);
    end
end

if arg.doTimeOnly % Make all vectors the correct length but zeros
    [rf, th, gx, gy, gz] = deal(zeros(nsamples,1));
else
    rf = rho.*exp(1i*th);
end

% plot
if arg.doDisplay
    T = (0:(numel(rho)-1))*raster/1000; % msec
    if ~arg.drawpause
        Tend = 1.01*T(find(any([rho th gx gy gz],2),1,'last')); %Find last non-zero value in any of the waveforms
    else
        Tend = 1.01*T(end);
    end
    
    gmax = arg.gmax; %5;  % Gauss/cm
    srho = arg.rhomax; %max(1.1*max(abs(rho(:))),0.05);
    lw = 1.5;
    subplot(511); plot(T, rho, 'LineWidth', lw); ylabel('|b1| (Gauss)'); axis([T(1) Tend -srho srho]);
    subplot(512); plot(T, th, 'LineWidth', lw);  ylabel('∠b1 (rad)'); axis([T(1) Tend -1.3*pi 1.3*pi]);
    subplot(513); plot(T, gx, 'LineWidth', lw);  ylabel('gx (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    %gmax = 1;  % Gauss/cm
    subplot(514); plot(T, gy, 'LineWidth', lw);  ylabel('gy (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    subplot(515); plot(T, gz, 'LineWidth', lw);  ylabel('gz (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    xlabel('msec');
end

return;

% EOF
