function [gx, gy, gz, gpre, esp, gx1, kz] = getcaipiepireadout(fov, imSize, Ry, Rz, Delta, gMax, slewRead, slewPre, raster, fbesp)
%
% Created 3D EPI CAIPI readout waveform.
% To generated single-echo GRE readout, use Ry = imSize(2)
%
% Inputs:
%  fov       [1 3]  cm
%  imSize    [1 3]  image volume matrix size
%  Ry        [1]    ky acceleration factor (can be non-integer)
%  Rz        [1]    kz acceleration factor (integer)
%  Delta     [1]    size of kz step (integer multiples of 1/fov(3))). Must be > 0.
%  gMax      Gauss/cm
%  slewRead  [1 3]  max slew along x/y/z gradient axis (Gauss/cm/ms)
%  slewPre   [1]
%  raster    ms
%  fbesp     [2] forbidden echo spacing range (ms)
%
% Outputs:
%  gx     x readout waveform (G/cm), without prephaser at beginning
%  gy     y readout waveform, without prephaser at beginning
%  gz     y readout waveform, without prephaser at beginning
%  gpre.x, gpre.y, gpre.z   x/y/z prephaser (to go to corner of kspace)
%  esp    echo spacing (ms)
%  gx1    waveform for one echo (G/cm)
%  kz     kz encoding indeces

% check inputs
if mod(Delta,1) | Delta < 1
    error('Delta must be positive integer');
end

nx = imSize(1); ny = imSize(2); nz = imSize(3);

dt = raster;
gamma = 4257.6;        % Hz/G

res = fov(1)/nx;       % spatial resolution (cm)
kmax = 1/(2*res);      % cycles/cm
area = 2*kmax/gamma;   % G/cm * sec (area of each readout trapezoid)
dkx = 1/fov(1);

dky = Ry/fov(2);       % ky spacing (cycles/cm)
dkz = 1/fov(3);        % kz spacing (cycles/cm) corresponding to Delta=1

etl = ny/Ry;

if mod(etl, 1)
    error('ny/Ry must be an integer')
end

% kz-encode indeces for one echo train
kz = [];
for ii = 1:Delta
    kz = [kz ii:Delta:Rz];
end
kz = repmat(kz, [1 ceil(etl/Rz)]);
kz = kz(1:etl);

% kz encoding blip amplitudes (multiples of dkz)
kzAmp = [diff(kz)];

% phase-encode blips
gyBlip = toppe.utils.trapwave2(dky/gamma, gMax, slewRead(2), dt);  % PE1
nyBlip = length(gyBlip);
gzBlip = toppe.utils.trapwave2(dkz*max(abs(kzAmp))/gamma, gMax, slewRead(3), dt);  % PE2
nzBlip = length(gzBlip);
nBlipMax = max([nyBlip nzBlip]);

% readout trapezoid (with ramp sampling)
mxg = min(1/(fov(1)*gamma*dt*1e-3), gMax);    % Gauss/cm
gx1 = toppe.utils.trapwave2(area, mxg, slewRead(1), dt);

% Reduce peak gradient until echo spacing is outside forbidden range
esp = length(gx1)*dt;  % echo spacing (ms)
if esp > fbesp(1) & esp < fbesp(2)
    for s = 1:-0.005:0.1
        mxg = s*mxg;
        gx1 = toppe.utils.trapwave2(area, mxg, slewRead(1), dt);
        if length(gx1)*dt > fbesp(2)
            esp = length(gx1)*dt;
            break;
        end
    end
end

% if needed, extend readout plateau to make room for turns
gx1orig = gx1;
nExtend = 0; % number of samples to add in middle of gx1
area1 = sum(gx1(ceil(nBlipMax/2):(end-ceil(nBlipMax/2))))*dt*1e-3;   % G/cm/s
while area1 < area
    nExtend = nExtend + 1;
    gx1 = [gx1orig(1:ceil(end/2)) max(gx1)*ones(1,nExtend) gx1orig((ceil(end/2)+1):end)];
    area1 = sum(gx1(ceil(nBlipMax/2):(end-ceil(nBlipMax/2))))*dt*1e-3;   % G/cm/s
end

% remove 0 at end in preparation for assembling into echo train
gx1 = gx1(1:(end-1));  

% echo train
gx = [];
for iecho = 1:etl
    gx = [gx gx1*(-1)^(iecho+1)];
end
gx = [gx(:); 0];  % waveform must end with 0

% add y and z blips
nt = length(gx1);  % number of samples in one echo
gy = 0*gx;
gz = 0*gx;
for ii = 1:(etl-1)
    iStart = ii*nt - nyBlip/2;
    iStop = iStart + nyBlip-1;
    gy(iStart:iStop) = gyBlip;

    iStart = ii*nt - nzBlip/2;
    iStop = iStart + nzBlip-1;
    gz(iStart:iStop) = gzBlip * kzAmp(ii)/max(abs(kzAmp));
end

% make length multiple of 4 
gx = toppe.makeGElength(gx(:));
gy = toppe.makeGElength(gy(:));
gz = toppe.makeGElength(gz(:));

% x/y/z prephasers.
% x prephaser a bit bigger since area of gx1 is a bit larger than 'area'
areax = sum(gx1)*dt*1e-3;   % G/cm * sec
gpre.x = toppe.utils.trapwave2(areax/2, gMax, slewPre, dt);
gpre.y = toppe.utils.trapwave2(area/2-dky/Ry/gamma/2, gMax, slewPre, dt);
gpre.y = [gpre.y(:); zeros(length(gpre.x)-length(gpre.y), 1)]; % make same length
%gpre = gpre(1:(end-1)); % remove 0 at end
gpre.x = toppe.makeGElength(gpre.x(:)); % make length multiple of 4
gpre.y = toppe.makeGElength(gpre.y(:));
gpre.z = gpre.y; % isotropic resolution

% plot k-space. Add prephaser for plotting purposes only.
gxf = [-gpre.x; gx];
gyf = [-gpre.y; gy];
gzf = [-gpre.y; gz];
kxp = gamma*dt*cumsum(gxf);
kyp = gamma*dt*cumsum(gyf);
kzp = gamma*dt*cumsum(gzf);
subplot(121); plot(kxp,kyp,'b.'); axis equal;
xlabel('kx (cycles/cm)'); ylabel('ky (cycles/cm)');
T = dt*1e3*(1:length(kxp));
subplot(122); hold off; plot(T,kxp,'r'); hold on; plot(T,kyp,'g'); plot(T,kzp,'b'); hold off;
legend('kx', 'ky', 'kz'); ylabel('cycles/cm'); xlabel('time (ms)');

% get time matrix

return

