function d = loadexchangedata(pfile, readout, din)
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

[ndat ncoils ndabslices ndabechoes ndabviews] = size(din);

nscans = ndabslices;
nz = ndabechoes;
nshots = ndabviews;

% get readout waveform
% see makeepi.m
[rf,gx,gy,gz,desc,paramsint16,paramsfloat,hdr] = toppe.readmod(readout);
nx = paramsint16(1); % decimation = 1 for this scan
ny = paramsint16(2);  % image matrix size is [nx ny nz]
nshots = paramsint16(3);
Ry = paramsint16(4);
nSampPerEcho = paramsint16(5) + paramsint16(6);
nramp = paramsint16(7);

etl = ny/nshots/Ry;  % echo train length

% sort data into [ny ny nz ncoils nscans] matrix
% size(din) = [ndat ncoils nscans nz nshots]
fprintf('Reshaping data matrix...');
d = zeros(nx, ncoils, nscans, nz, ny);
for ish = 1:nshots
    for e = 1:etl  % loop over echoes
        iy = ish + (e-1)*nshots;
        iStart = (e-1)*nSampPerEcho + nramp + 1;
        iStop = iStart + nx - 1;
        d(:, :, :, :, iy) = din(iStart:iStop, :, :, :, ish);   
    end
end
d = permute(d, [1 5 4 2 3]);
fprintf('done\n');

return

% slice to write to .mat file
sl = 5;   % middle slice

% recon
fprintf('recon..');
% recon complex coil SPGR images (for B1+ and B0 estimation)

for s = 1:nscans
    fprintf('.');
    [spgr(:,:,:,:,s) tmp] = toppe.utils.ift3(d(:,:,:,:,s));
end
fprintf('\t');
spgr = spgr((end/2-nx/2):(end/2+nx/2-1), :, sl, :, :);
spgr = squeeze(spgr);  % [nx ny nz ncoils 2]

if scanType == 1
    save('bs.mat', '-v7.3', 'spgr');
    return;
else
    save(sprintf('spgr,design%d.mat', scanType), '-v7.3', 'spgr');
end

% recon root-sum-of-squares STFR images
clear stfr
for s = 3:size(d,5)
    fprintf('.');
    [tmp stfr(:,:,:,s-2)] = toppe.utils.ift3(d(:,:,:,:,s));
end
stfr = stfr((end/2-nx/2):(end/2+nx/2-1), :, sl, :);
stfr = squeeze(stfr);  % [nx ny nz nscans-2]
im(stfr, [0 0.3*max(stfr(:))]);
fprintf(' done\n');

save(sprintf('stfr,design%d.mat', scanType), '-v7.3', 'stfr');

return
