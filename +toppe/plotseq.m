function [rf, gx, gy, gz, rf1, gx1, gy1, gz1, tdelay] = plotseq(nstart, nstop, system, varargin)
% function [rf, gx, gy, gz, rf1, gx1, gy1, gz1, tdelay] = plotseq(nstart, nstop, system, varargin)
%
% Display pulse sequence, as specified in modules.txt, scanloop.txt, and timing.txt.
% Sequence timing is exact (at least that's the intent), i.e., it
% tries to match the timing on the scanner exactly.
% See resources/timing.svg for a detailed sequence timing diagram
%
% Inputs:
%   nstart,nstop       first and last startseq calls (as specified in scanloop.txt)
%   system             struct specifying hardware system info, see systemspecs.m
%
% Input options:
%   loopFile           default: 'scanloop.txt'
%   loopArr            scan loop array (see readloop.m). Default: read from loopFile.
%   moduleListFile     default: 'modules.txt'
%   mods               Structure containing .mod file contents (see ./utils/readModulesListFile.m).
%                      Use in playseq.m to speed up display.
%                      Default: get values from moduleListFile.
%   doDisplay          true (default) or false
%   doTimeOnly         Returns outputs as zeros (but with correct length)
%                      to speed up calculation of scan time.
%                      False (default) or true
%   drawpause          (boolean) include pauses (textra) or not
%   gmax               display limit, Gauss/cm
%   rhomax             display limit, Gauss
%   printTime          False (default) or true
%
% Outputs:
%   rf               Complex RF waveform (Gauss)
%   gx,gy,gz         Gauss/cm
%   rf1/gx1/...      Values from most recent startseq() call (used in, e.g., ge2seq.m and playseq.m)
%   tdelay           Delay after end of module waveform.
%                    Determined by duration in modules.txt AND by textra (column 14) in scanloop.txt
%                    Used in ge2seq.m to set delay block, and in playseq.m.

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
    cores = toppe.tryread(@toppe.readmodulelistfile, arg.moduleListFile);
else
    cores = arg.mods;
end

%% Initialize counter and turn off display if we're only doing timings
if arg.doTimeOnly
    nsamples = 0;
    arg.doDisplay = false;
    for ic = 1:size(cores,2) % Compute table of core durations
        core_size(ic) = size(cores{ic}.gx(:,1),1);
    end
end

%% build sequence. each sample is 4us.
rho = []; th = []; gx = []; gy = []; gz = [];
dt = 4;  % us
max_pg_iamp = 2^15-2;
for it = nstart:nstop
    ic = loopArr(it,1);   % core id
    if cores{ic}.hasRF
        ia_rf = loopArr(it,2);
    else
        ia_rf = 0;
    end
    ia_th = loopArr(it,3);
    ia_gx = loopArr(it,4);
    ia_gy = loopArr(it,5);
    ia_gz = loopArr(it,6);
    
    % start of core
    if cores{ic}.hasRF
        start_core = max(system.start_core_rf - dt*cores{ic}.npre, 0);
    elseif cores{ic}.hasDAQ
        start_core = max(system.start_core_daq - dt*cores{ic}.npre, 0);
    else
        start_core = system.start_core_grad;
    end

    % number of discarded samples at end of RF/ADC window
    nChopEnd = cores{ic}.res - cores{ic}.rfres - cores{ic}.npre;

    % delay past gradient due to RF/DAQ window
    if cores{ic}.hasRF
        coredel = max(system.myrfdel - dt*nChopEnd, 0);
    elseif cores{ic}.hasDAQ
        coredel = max(system.daqdel - dt*nChopEnd, 0);
    else
        coredel = 0;
    end

    % Mimimum core duration (us).
    % Should be identical to the CV 'mindur' in the EPIC code.
    mindur = start_core + cores{ic}.wavdur + system.timetrwait + coredel + system.tminwait + system.timessi;

    % silence at end of core
    tdelay = max(cores{ic}.dur - mindur, 0); 
    if size(loopArr,2)>13
        tdelay = tdelay + loopArr(it,14);  % add textra
    end
    
    waveform = loopArr(it,16); % waveform index
    
    if arg.doTimeOnly % Calculate the length of one waveform and add it to our sample counter
        gxlength = round((start_core)/dt) + core_size(ic) + round((system.timetrwait+system.timessi+coredel)/dt);
        nsamples = nsamples + gxlength + round(tdelay/dt);
    else % Calculate RF and gradients as normal
        % get gradients
        gxit = ia_gx/max_pg_iamp*cores{ic}.gx(:,waveform);
        gyit = ia_gy/max_pg_iamp*cores{ic}.gy(:,waveform);
        gzit = ia_gz/max_pg_iamp*cores{ic}.gz(:,waveform);

        % apply 3d rotation matrix 
        % (which also accounts for any in-plane 2D rotation, i.e., 'phi' in write2loop.m)
        Rv = loopArr(it,17:25)/max_pg_iamp;  % stored in row-major order
        R = reshape(Rv, 3, 3);
        G = R * [gxit(:)'; gyit(:)'; gzit(:)'];
        gxit = G(1,:)';
        gyit = G(2,:)';
        gzit = G(3,:)';
        
        % build waveforms for this startseq call
        rho1 = [zeros(round((start_core+coredel)/dt),1); ...
                ia_rf/max_pg_iamp*  abs(cores{ic}.rf(:,waveform)); ...
                zeros(round((system.tminwait+system.timetrwait+system.timessi)/dt),1)];
        th1  = [zeros(round((start_core+coredel)/dt),1); ...
                ia_th/max_pg_iamp*angle(cores{ic}.rf(:,waveform)); ...
                zeros(round((system.tminwait+system.timetrwait+system.timessi)/dt),1)];
        gx1  = [zeros(round((start_core)/dt),1); ...
                gxit(:); ...
                zeros(round((system.tminwait+system.timetrwait+system.timessi+coredel)/dt),1)];
        gy1  = [zeros(round((start_core)/dt),1); ...
                gyit(:); ...
                zeros(round((system.tminwait+system.timetrwait+system.timessi+coredel)/dt),1)];
        gz1  = [zeros(round((start_core)/dt),1); ...
                gzit(:); ...
                zeros(round((system.tminwait+system.timetrwait+system.timessi+coredel)/dt),1)];
        
        % apply RF phase offset
        if cores{ic}.hasRF
            th1 = th1 + loopArr(it,12)/max_pg_iamp*pi;
            th1 = angle(exp(1i*th1));   % wrap to [-pi pi] range
        end

        % add to running waveform
        rho = [rho; rho1; zeros(round(tdelay/dt),1)];
        th  = [th;  th1;  zeros(round(tdelay/dt),1)];
        gx  = [gx;  gx1;  zeros(round(tdelay/dt),1)];
        gy  = [gy;  gy1;  zeros(round(tdelay/dt),1)];
        gz  = [gz;  gz1;  zeros(round(tdelay/dt),1)];
    end
    
    if arg.printTime
        fprintf(1, 'it %d: mindur = %d us, rf t = %d us, grad t = %d us\n', it, mindur, numel(rho)*dt, numel(gx)*dt);
    end
end

if arg.doTimeOnly % Make all vectors the correct length but zeros
    [rf, th, gx, gy, gz] = deal(zeros(nsamples,1));
else
    rf = rho.*exp(1i*th);
    rf1 = rho1.*exp(1i*th1);     % waveforms in last module, without the delay after it (if any)
end

if length(rf) ~= length(gx)
keyboard
    msg = ['RF and gradient waveform durations are not equal length. ' ... 
        'This is likely due to system.rfdel or system.daqdel not being ' ...
        'on a 4us boundary.'];
    error(msg);
end

% plot
if arg.doDisplay
    T = (0:(numel(rho)-1))*dt/1000; % msec
    if ~arg.drawpause
        Tend = 1.01*T(find(any([rho th gx gy gz],2),1,'last')); %Find last non-zero value in any of the waveforms
    else
        Tend = 1.01*T(end);
    end
    
    gmax = arg.gmax; %5;  % Gauss/cm
    srho = arg.rhomax; %max(1.1*max(abs(rho(:))),0.05);
    lw = 1.5;
    subplot(511); plot(T, rho, 'LineWidth', lw); ylabel('rho (Gauss)'); axis([T(1) Tend -srho srho]);
    subplot(512); plot(T, th, 'LineWidth', lw);  ylabel('theta (rad)'); axis([T(1) Tend -1.3*pi 1.3*pi]);
    subplot(513); plot(T, gx, 'LineWidth', lw);  ylabel('gx (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    %gmax = 1;  % Gauss/cm
    subplot(514); plot(T, gy, 'LineWidth', lw);  ylabel('gy (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    subplot(515); plot(T, gz, 'LineWidth', lw);  ylabel('gz (G/cm)'); axis([T(1) Tend -1.05*gmax 1.05*gmax]);
    xlabel('msec');
end

return;

% EOF
