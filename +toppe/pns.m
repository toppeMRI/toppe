function [PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,coil,varargin)
% function [PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,chronaxie,rheobase,alpha,gdt,plt,varargin)
%PNS Convolution model to calculate peripheral nerve stimulation
% [PThresh,pt,PTmax,gmax,smax,t,f] = pns(grad,chronaxie,rheobase,alpha,gdt,plt)
%                                              CV par on scanner
%   PThresh  Threshold of PNS          [%]     (cfdbdtper)
%        pt  Threshold of different axis (uncombined)
%     PTmax  Maximum PNS threshold (combined)
%      gmax  Maximum gradient strength (combined)
%      smax  Maximum slewrate (combined)
%
% Inputs:
%      grad  gradient matrix (dimensions,points,repetitions) [T/m]
%            (for 2D also dim=1 + complex data possible)
%            also possible: (points,repetitions,dimensions)
%      coil  'xrm', 'xrmw', 'whole', 'zoom' (take values from list below)
%
% Input options:
%      print    (bool) print a few relevant output values to console
%      plt      (bool) plot
%      gdt      (sec) Gradient update time (default=4d-6)

% Old inputs:
% chronaxie  Chronaxie time constant   [s]     (cfrfact [us])
%  rheobase  Rheobase scale factor     [T/s]   (cfrinf [T/s])
%     alpha  Effective coil length     [m]     (cfdbdtdx/y/z [cm])
%            SRmin=rheobase/alpha
%       gdt  Gradient update time      [s] (default=4d-6)
%       plt  Plotting + print output (default=true)
%
% Scanner  Gradient coil   chronaxie rheobase alpha  gmax  smax
% MR750w   XRMW            360d-6    20.0     0.324  33    120
% MR750    XRM             334d-6    23.4     0.333  50    200
% HDx      TRM WHOLE       370d-6    23.7     0.344  23    77
% HDx      TRM ZOOM        354d-6    29.1     0.309  40    150
% UHP      HRMB            359d-6    26.5     0.370  100   200
%
% values on scanner from /w/config/Scandbdt.cfg or GRSubsystemHWO.xml
% (e.g., /export/home/mx/host/config/current/GRSubsystemHWO.xml)
% (alpha = EffectivedBdTlength<X,Y,Z>/100)
%
% Literature: "Introduction to dB/dt" presentation by Toni Linz; 6/8/06
%             DIN EN 60601-2-33
% 2/2014 Rolf Schulte
%
% Modified by Jon-Fredrik Nielsen for TOPPE 2018-2019
% 2-Oct-2021  JFN  Added PNS model parameters for HRMB (UHP) system

import toppe.utils.*       % vararg_pair

arg.print = true;
arg.gdt = 4d-6;       % sec
arg.plt = true;
arg = vararg_pair(arg, varargin);

if nargin<1, help(mfilename); return; end

model = 1;         % 1=convolution in time domain; 2=multiplication in FD

%% input checking
if length(size(grad))==3,
    [n1,n2,n3] = size(grad);
    if ((n1>n2) && (n1>n3))
        grad = shiftdim(grad,2);   % JFN: (dimensions,points,repetitions)
    end
end
if ~isreal(grad),
    if size(grad,1)~=1, error('grad complex and size(grad,1)~=1'); end
    grad = [real(grad) ; imag(grad)];
end
if any(isnan(grad(:))), warning('pns:nan','grad contains NaN'); end
if any(isinf(grad(:))), warning('pns:inf','grad contains inf'); end

switch lower(coil)
    case 'xrmw',  chronaxie=360d-6; rheobase=20.0; alpha=0.324;
    case 'xrm',   chronaxie=334d-6; rheobase=23.4; alpha=0.333;
    case 'whole', chronaxie=370d-6; rheobase=23.7; alpha=0.344;
    case 'zoom',  chronaxie=354d-6; rheobase=29.1; alpha=0.309;
    case 'hrmb',  chronaxie=359d-6; rheobase=26.5; alpha=0.370;
    otherwise, error('gradient coil (%s) unkown',chronaxie);
end

if ~exist('rheobase','var'), error('Provide rheobase'); end
if ~exist('alpha','var'), error('Provide alpha'); end
if ~exist('gdt','var'), gdt = []; end
if isempty(gdt), gdt = 4d-6; end
if (chronaxie>600d-6) || (chronaxie<200d-6),
    warning('pns:chronaxie','chronaxie=%g; typical values in help',chronaxie);
end
if (rheobase>30) || (rheobase<20),
    warning('pns:rheobase','rheobase=%g; typical values in help',rheobase);
end
if (alpha>0.4) || (alpha<0.3),
    warning('pns:alpha','alpha=%g; typical values in help',alpha);
end

%% zerofilling for FFT (must decay to zero to avoid wiggles in front)
if model==2,
    if isodd(size(grad,2)), grad(:,size(grad,2)+1,:) = 0; end
    zf2 = ceil(7d-3/gdt/2)*2;
    dozf = false;
    if size(grad,2)<zf2,
        dozf = true;
    else
       if any(any(grad(:,end-zf2+1:end,1)~=0)), dozf = true; end
    end
    if dozf, grad(:,size(grad,2)+zf2,:) = 0; end
end

if length(size(grad))>3, error('length(size(grad))>3'); end
if length(size(grad))<2, error('length(size(grad))<2'); end
[n1,n2,n3] = size(grad);
if n1>3, error('size(grad,1)>3'); end
if n2==1, error('size(grad,2)==1'); end
if n3>n2, warning('pns:n3_n2','size(grad,3)>size(grad,2)'); end


%% misc parameters
SRmin = rheobase/alpha;      % min slewrate causing stimulation; Slide 3 Eq 3+4 [T/m/s]
SR = [zeros(n1,1,n3),diff(grad,[],2)]/gdt;  % instantaneous slewrate [T/m/s]

switch model
    case 1,  % convolution in time domain (identical to calcpns)
        pt = zeros(n1,n2,n3);
        for l3=1:n3,
            for l2=1:n2,
                % Slide 34 Eq 30
                pt(:,l2,l3) = 100*gdt*chronaxie/SRmin*...
                    sum(SR(:,(1:l2),l3)./(chronaxie+repmat((l2:-1:1)-0.5,[n1,1])*gdt).^2,2);
            end
        end
    case 2,  % multiplication in freq domain
        % small difference exist because time is shifted 1/2*gdt in time

        decay = fftshift(fft(1./(((0:n2-1)+0.5)*gdt+chronaxie).^2));
        % decay = fftshift(fft(1./(((0:n(2)-1))*gdt+chronaxie).^2));
        SR_fd = fftshift(fft(SR,[],2),2);
        pt = 100*gdt*chronaxie/SRmin*...
            ifft(ifftshift(SR_fd.*repmat(decay,[n1,1,n3]),2),[],2);
    otherwise, error('model (%g) unknown',model);
end

if any(isnan(pt(:))), warning('pns:nan','pt contains NaN'); end
if any(isinf(pt(:))), warning('pns:inf','pt contains inf'); end

PThresh = sqrt(sum(pt.^2,1));      % percent threshold; Slide 35; Eq 31


%% plotting
if arg.plt
    t = (0:n2-1)*gdt*1d3;            % time [ms]
    %figure; %(13); 
    clf;
    rsos_grad = sqrt(sum(grad(:,:,1).^2,1));
    subplot(2,2,1); plot(t,grad(:,:,1)*1d3,'',t,rsos_grad*1d3,'r--');
    xlabel('time [ms]'); ylabel('grad [mT/m]');
    tmp = ceil(max(abs(rsos_grad))*1050);
    %grid on; axis([0 t(end) -tmp tmp]);
    
    rsos_SR = sqrt(sum(SR(:,:,1).^2,1));
    subplot(2,2,3); plot(t,SR(:,:,1),'',t,rsos_SR,'r--');
    xlabel('time [ms]'); ylabel('slewrate [T/m/s]');
    tmp = ceil(max(abs(rsos_SR))*1.05);
    %grid on; axis([0 t(end) -tmp tmp]);
    
    subplot(2,2,4); 
    plot(t,pt(:,:,1),'',t,PThresh(:,:,1),'r--',...
        [t(1) t(end)],[100 100],'m:',[t(1) t(end)],[80 80],'m:',...
        [t(1) t(end)],-[100 100],'m:',[t(1) t(end)],-[80 80],'m:');
    xlabel('time [ms]'); ylabel('PNS threshold [%]');
    tmp = 1.05*max([PThresh(:,:,1) 100]);
    grid on; axis([0 t(end) -tmp tmp]);
    
    if model==2,
        f = (-n2/2:(n2-1)/2)/n2/gdt*1d-3;
        subplot(2,2,2);
        plot(f,abs(SR_fd(:,:,1))/max(max(abs(SR_fd(:,:,1)))),'',...
            f,abs(decay)/max(abs(decay(:))),'r');
        xlabel('freq [kHz]'); ylabel('FFT of decay + slewrate');
        grid on; axis([-5 5 -0.05 1.05]);
			legend
    end
end

%% print output
if arg.print,
    fprintf('\nchronaxie = %g [us]\n',chronaxie*1d6);
    fprintf('rheobase  = %g [T/s]\n',rheobase);
    fprintf('effective coil length: alpha = %g [cm]\n',alpha*1d2);
    fprintf('Stimulation slewrate:  SRmin = %g [T/m/s]\n',SRmin);
    fprintf('Gradient update time:    gdt = %g [us]\n',gdt*1d6);
    fprintf('Maximum gradient strength    = %g [mT/m]\n',max(grad(:))*1d3);
    fprintf('Maximum slewrate             = %g [T/m/s]\n',max(SR(:)));

    if n3==1, fprintf('Maximum PThresh = %.4g [%%]\n',max(abs(PThresh)));
    else
        fprintf('Maximum PThresh = \n');
        for l3=1:n3, fprintf('\tl3=%g: %.4g [%%]\n',l3,max(abs(PThresh(:,:,l3)))); end
        fprintf('\tall: %.4g [%%]\n',max(abs(PThresh(:))));
    end
    
    if max(PThresh)>80,
        if max(PThresh)>100,
            fprintf('Warning: PThresh exceeding first controlled mode (100%%)!!!\n');
        else
            fprintf('Warning: PThresh exceeding normal mode (80%%)!\n');
        end
    end
end

%% misc other output arguments
if nargout>2, PTmax = max(PThresh); end
if nargout>3, 
    if ~exist('rsos_grad','var'), rsos_grad = sqrt(sum(grad(:,:,1).^2,1)); end
    gmax = max(rsos_grad); 
end
if nargout>4, 
    if ~exist('rsos_SR','var'), rsos_SR = sqrt(sum(SR(:,:,1).^2,1)); end
    smax = max(rsos_SR); 
end

return;


function val = isodd(n)

val = mod(n,2) > 0

return
