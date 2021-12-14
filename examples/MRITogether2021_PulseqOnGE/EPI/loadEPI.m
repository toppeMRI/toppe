function [d, ims] = loadEPI(pfile, readout, din)
% Run this from the location of the Pfile, e.g.,
% /mnt/storage/jfnielse/exchange/...

if nargin < 3
% load full data matrix
din = toppe.utils.loadpfile(pfile);  % [ndat ncoils nslices nechoes nviews]
din = flipdim(din, 1); % !! Data seems to always be flipped along readout (do check though)

% permute
% scans store slice in 'echo' index, since max echoes = 16
%din = permute(d, [1 5 4 2 3]);  % [ndat nshots nz ncoils nscans]
end

[ndat ncoils nz nscans nshots] = size(din);

% get readout waveform
% see makeepi.m
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readout);
nx = paramsint16(1); % decimation = 1 for this scan
ny = paramsint16(2);  % image matrix size is [nx ny nz]
nshots = paramsint16(3);
Ry = paramsint16(4);
nSampPerEcho = paramsint16(5) + paramsint16(6);
nramp = paramsint16(7);

% echo train length
etl = ny/nshots/Ry;  

% sort data into [nx ny nz ncoils nscans] matrix
% size(din) = [ndat ncoils nscans nz nshots]
fprintf('Reshaping data matrix...');
d = zeros(nx, ncoils, nz, nscans, ny);
for ish = 1:nshots
    for e = 1:etl  % loop over EPI echoes
        iy = ish + (e-1)*nshots;
        iStart = (e-1)*nSampPerEcho + nramp + 1;
        iStop = iStart + nx - 1;
        d(:, :, :, :, iy) = din(iStart:iStop, :, :, :, ish);   
    end
end
d = permute(d, [1 5 3 2 4]);   % [nx ny nz ncoils nscans]
fprintf(' done\n');

% Reconstruct
ims = zeros(nx, ny, nz, nscans);   % root-sum-of-squares coil-combined images
for isc = 1:nscans
    [~, ims(:,:,:,isc)] = toppe.utils.ift3(d(:,:,:,:,isc));
end
figure
im(ims(:,:,:,1));
title('3D image volume, scan 1 (flip angle = 5 deg)')
input('Press enter to proceed.');

% Extract one slice and display images for all flip angles
figure
sl = 6;  % slice
ims = squeeze(ims(:,:,sl,:));
im(ims)
title('Slice 6, flips = 5:5:40')
input('Press enter to proceed.');

% Plot signal vs flip angle in an ROI
figure
imagesc(ims(:,:,1));
title('Draw ROI');
roi = roipoly;
for isc = 1:nscans
    tmp = ims(:,:,isc);
    s(isc) = mean(tmp(roi));
end
flip = 5:5:40; % degrees
plot(flip, s, 'bo');
title('Signal vs flip angle');
xlabel('flip angle (degrees)');
ylabel('SPGR signal (a.u.)');
hold on;

% Add the theoretical SPGR-vs-flip angle curve 
TR = 26.5; % ms
T1 = 1000;  % FBIRN phantom (guess)
E1 = exp(-TR./T1);
flip2 = 1:40;
s_theory = (1-E1).*sind(flip2) ./ (1-cosd(flip2).*E1);
hold on;
plot(flip2, s_theory*max(s)/max(s_theory), 'r-');
legend('observed', sprintf('SPGR/FLASH signal equation (T1 = %d ms)', T1));

ernst = acos(E1)/pi*180;

return

