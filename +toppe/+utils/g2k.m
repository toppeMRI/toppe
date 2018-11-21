function [kx,ky] = g2k(g,nint,OFF)
% function [kx,ky] = g2k(g,nint,OFF)
%
% INPUTS:
%    g:     [N 2] array containing Gx and Gy in Gauss/cm
%    nint:  number of leaf rotations (default=1)
%    OFF:   [1 2] sub-sample shift in gx and gy (% of a sample). Default: [50 50]

g = complex(g(:,1),g(:,2));

if ~exist('nint','var')
	nint=1;
end

% calculate k-space from gradients
gambar = 4.258e3;   % Hz/gauss 
gts = 4e-6;         % gradient sample duration (dwell time) (sec)
if ~exist('OFF','var')
	OFF = [50 50];
end
R = 100;
gtemp = interp(g,R);
gamma = 4.2576;      % kHz/Gauss
dt = 4e-3;           % msec

ktemp = gamma*dt*cumsum(gtemp)/R;   % cycles/cm
%ktemp = gamma*dt*complex(cumtrapz(real(gtemp)), cumtrapz(imag(gtemp)))/R;
ktemp = [real(ktemp) imag(ktemp)];

k = [ktemp(OFF(1):R:size(ktemp,1),1) ktemp(OFF(2):R:size(ktemp,1),2)];
k = complex(k(:,1),k(:,2));
%K = cumsum([ g])*gts*gambar;      % cycles/cm
%Ky = cumsum([ g(:,2)])*gts*gambar;

kspace(:,1) = real(k);
kspace(:,2) = imag(k);

kxi = kspace(:,1);
kyi = kspace(:,2);
Gx = real(g);
Gy = imag(g);

tsamp = gts;
kxt=interp1([0:4e-6:4e-6*length(kxi)-4e-6],kxi,[0:tsamp:4e-6*length(kxi)-tsamp])';
kyt=interp1([0:4e-6:4e-6*length(kyi)-4e-6],kyi,[0:tsamp:4e-6*length(kyi)-tsamp])';

gxt=interp1([0:4e-6:4e-6*length(Gx)-4e-6],Gx,[0:tsamp:4e-6*length(Gx)-tsamp])';
gyt=interp1([0:4e-6:4e-6*length(Gx)-4e-6],Gy,[0:tsamp:4e-6*length(Gx)-tsamp])';

nk = length(kxt)-2;   

kx = zeros(nk,nint);
ky = zeros(nk,nint);
kxo = zeros(nk,1);
kyo = zeros(nk,1);

kxo = kxt(1:nk);
kyo = kyt(1:nk);

%rotate matrix (in-plane) for proper orientation
rotamount = 0;
phir = -rotamount*pi/2;
kxop = kxo*cos(phir) - kyo*sin(phir);
kyop = kyo*cos(phir) + kxo*sin(phir);

%fprintf('%s.m: Performing %d rotations\n',mfilename,nint);
kx(:,1) = kxop;
ky(:,1) = kyop;
phi = 2*pi/nint;
for ii = 1:(nint-1)
     kx(:,ii+1) = kxop*cos(ii*phi) - kyop*sin(ii*phi);
     ky(:,ii+1) = kyop*cos(ii*phi) + kxop*sin(ii*phi);
end

return;
