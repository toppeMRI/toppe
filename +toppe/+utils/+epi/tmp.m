function [gx,gy,gz] = makeepi(fov, nPix, nEcho, varargin)
% function [gx,gy,gz] = makeepi(fov, nPix, nEcho, varargin)
%
% Make non-flyback multi-echo readout. Play out with TOPPE pulse sequence.
%
% INPUTS:
%  fov     - cm
%  npix    - number of pixels 
% Options:
%  necho   - number of echoes
%
% $Id: tmp.m,v 1.1 2018/10/28 18:03:35 jfnielse Exp $
% $Souce: $

isbalanced = 0;

dt = 4e-3;               % gradient/daq sample duration (msec)
mxg = 4.9;
mxs = 10;              % go easy on gradients


%% x gradient

% set readout amplitude and number of datapoints to acquire, based on oprbw
gamma = 4.2575;                     % kHz/Gauss
oprbw = 125/4;
scalefactor = 125/oprbw;
g = (1/dt)/(gamma*fov)/scalefactor;             % Gauss/cm
npix = npix*scalefactor;
if g > mxg
	g = mxg;
	fprintf(1, 'max g reduced to %.2f G/cm \n', g);
end

% readout trapezoid
gxro = g*ones(1,npix);                         % plateau of readout trapezoid
s = mxs*dt;
ramp = s:s:(g-s) ;  
%gxro = [ramp g*ones(1,npix)];
gxro = [ramp g*ones(1,npix) fliplr(ramp)];

% x prewinder. make sure res_kpre is even. 
area = sum(gxro)*4e-6;
gxprew = -trapwave(area/2, 4e-6, mxg, mxs*1e3);    % 1 mm resolution phase-encode blip
gxprew = makeevenlength(gxprew);

% balancing gradient between echoes (x)
gbal = -trapwave(area, 4e-6, mxg, mxs/sqrt(2)*1e3);    % sqrt(2) since playing both x and z gradients

gx = [gxprew gxro];
for ii = 1:(necho-1)
	gx = [gx gbal gxro];
end
if isbalanced
	gx = [gx gxprew];
end
gy = 0*gx;
gz = 0*gx;

T = linspace(dt/2,length(gx)*dt-dt/2,length(gx));
figure; plot(T, gx, 'b');
hold on; plot(T, gy, 'r'); plot(T, gz, 'g');
xlabel('time (msec)');

% write to readout.wav
if necho > 1
	fname = sprintf('multiecho-fov%d-npix%d-necho%d.wav',fov,npix/scalefactor,necho);
else
	%fname = sprintf('fov%d-npix%d-isbalanced%d.wav',fov,npix,isbalanced);
	fname = sprintf('fov%d-npix%d-oprbw%.3f.wav',fov,npix/scalefactor,oprbw);
end
npre = 100;
rfres = npix;
mat2wav(0.01*ones(size(gx(:))),0*gx(:),gx(:),gy(:),gz(:),90,fname,'multi-echo readout waveform',npre,rfres);
plotwav(fname);

% done
 
return;

% EOF
