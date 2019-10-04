function plotrfspectrum(rf, fmt)
%
% Inputs
%  rf       [ntp 1]      assumes 4us raster (sample) time
%  fmt      plot format string, e.g., 'bo'

zpad = 500;
ft = fftshift(fft(fftshift([zeros(zpad,1); rf; zeros(zpad,1)])));
dt = 4e-6;             % sec
fov = length(ft)*dt;   % sec
df = 1/fov/1000;       % kHz
ft = ft((end/2-50):(end/2+50));
fmax = (length(ft)/2*df - df/2);
F = linspace(-fmax,fmax,length(ft));
plot(F, abs(ft), fmt)
xlabel('kHz');
ylabel('a.u.');

return;
