function [g,seq] = makesosreadout(sys, N, FOV, nLeafs, varargin)
% Make stack-of-spirals readout.mod file (fully sampled)
%
% function g = makesosreadout(sys, N, FOV, nLeafs, varargin)
%
% Inputs:
%   FOV        [x y z] field-of-view (cm)
%   N     [nx ny nz] reconstructed image matrix size
%   nLeafs     number of spiral leafs 
% Options:
%   dsamp      Number of samples for fully-sampled center (when spiralDesign = 'genspivd2'). Default: 500.
%   Router     Undersampling factor outside dense center (when spiralDesign = 'genspivd2'). Default: 2.
%   spiralDesign    'genspivd2' (default) or 'vds'
%   maxGrad    Max gradient amplitude. Default: system.maxGrad
%   inout      'in' or 'out' (default) for spiral-in/out
%   ofname     output file name. Default: 'readout.mod'
%   rewDerate  derate slew during rewinder and z phase-encode blip by this factor, to reduce PNS. See pns.m. Default: 0.8.
% Outputs:
%   g          [nt 3]       [gx gy gz] gradients (G/cm)
%   seq        struct       misc design parameters that may be useful (e.g., seq.sampWin)
%

import toppe.*
import toppe.utils.*
import toppe.utils.spiral.*
import toppe.utils.spiral.mintgrad.*

%% parse and check inputs
% Defaults
arg.dsamp = 500;
arg.Router = 2;
arg.spiralDesign = 'genspivd2';  % use either genspivd2.m or vds.m to design spiral
arg.maxGrad = sys.maxGrad;
arg.inout   = 'out';
arg.ofname  = 'readout.mod';
arg.rewDerate = 0.8;

%arg = toppe.utils.vararg_pair(arg, varargin);
arg = vararg_pair(arg, varargin);

if FOV(1) ~= FOV(2)
	error('Anisotropic x/y FOV not supported.');
end
if N(1) ~= N(2)
	error('Anisotropic x/y matrix not supported.');
end

%% struct to be returned
seq.FOV = FOV;
seq.N = N;
seq.nLeafs = nLeafs;
seq.system = sys;

maxSlew = 0.999*sys.maxSlew;

%% design spiral waveform (balanced)
npix = N(1);

doreverse= 0;
dovardens = 0;

rmax = npix/(2*FOV(1));   % max k-space radius

% genspivd2 returns complex g
if strcmp(arg.spiralDesign, 'genspivd2')
    FOVvd = FOV(1)/nLeafs;
    xresvd = N(1)/nLeafs;
    mxg = sys.maxGrad;   % G/cm
    mxslew = 10*sys.maxSlew;   % T/m/s
    [g] = toppe.utils.spiral.genspivd2(FOVvd, xresvd, arg.Router, mxg, mxslew, arg.dsamp);
else
    [g] = vds(maxSlew*1e3, sys.maxGrad, sys.raster, nLeafs, FOV(1), 0, 0, rmax);
    %cmd = sprintf('Created with Brian Hargreaves'' code: [k,g] = vds(%d,%d,4e-6,%d,fov,0,0,rmax);', maxSlew, sys.maxGrad, nLeafs);
end

g = [0; 0; g(:)];  % add a couple of zeroes to make sure k=0 is sampled
nsamp = length(g);

% make balanced
gx = makebalanced(real(g(:)), 'maxSlew', arg.rewDerate*maxSlew/sqrt(2));  
gy = makebalanced(imag(g(:)), 'maxSlew', arg.rewDerate*maxSlew/sqrt(2));   

% make same length
n = max(length(gx), length(gy));
gx = [gx; zeros(n-length(gx), 1)];
gy = [gy; zeros(n-length(gy), 1)];

% partition (kz) encoding trapezoid
gzamp = (1/sys.raster)/(sys.gamma*FOV(3));     % Gauss/cm
zarea = gzamp*N(3)*sys.raster;                   % Gauss/cm*sec
gpe = -trapwave2(zarea/2, sys.maxGrad, arg.rewDerate*maxSlew, sys.raster*1e3);

% put kz trapezoid and spiral together
gxtmp = gx;
gx = [0*gpe(:); zeros(2,1);   gx(:);    0*gpe(:)];
gy = [0*gpe(:); zeros(2,1);   gy(:);    0*gpe(:)];
gz = [  gpe(:); zeros(2,1); 0*gxtmp(:);  -gpe(:)];

seq.sampWin = (length(gpe)+5):(length(gpe)+5+nsamp-1);

% spiral in or out 
if strcmp(arg.inout, 'in')
	gx = flipud(gx);
	gy = flipud(gy);
	seq.sampWin = fliplr(length(gx)-seq.sampWin);
end

% make sure duration is on 4-sample (16us) boundary
gx = toppe.makeGElength(gx);
gy = toppe.makeGElength(gy);
gz = toppe.makeGElength(gz);

% write to .mod file
writemod(sys, 'gx', gx, 'gy', gy, 'gz', gz, 'ofname', arg.ofname, 'desc', 'stack-of-spirals readout module');

% return gradients for one leaf (to be rotated and kz-blipped in scanloop.txt)
g = [gx gy gz];

return;







if dovardens
	type = 'vds';
else
	type = 'unif';
end
if doreverse
	g = flipdim(g,1);  % reverse spiral 
	fname = sprintf('g-reverse-%s-nl%d-fov%d-npix%d-%s.mod', type, nLeafs, fov, npix, date);
else
	fname = sprintf('g-%s-nl%d-fov%d-npix%d.mod', type, nLeafs, fov(1), npix);
end
fprintf(1,'gradient duration is %.2f ms \n', 4e-3*length(g));
desc = sprintf('spiral readout .mod file for toppe\n%s\n',cmd);
%toppe.writemod('gx', gx, 'gy', gy, 'ofname', fname, 'desc', desc);

%Combines spiral readout and z phase-encode blips.

% readout for second DESS echo (spiral-in)
gx2 = flipud(gx1);
gy2 = flipud(gy1);
gz2 = flipud(gz1);

% write to readout.mod and plot
gx = [gx1, gx2];
gy = [gy1, gy2];
gz = [gz1, gz2];
gx = makeGElength(gx);   % make divisible by 4
gy = makeGElength(gy);
gz = makeGElength(gz);
rhfrsize = size(gx, 1);
%mat2mod(0.01*ones(size(gx)),0*gx,gx,gy,gz,90,'readout.mod','1mm isotropic DESS stack-of-spirals readout',0, [0 rhfrsize], false);
writemod('gx', gx, 'gy', gy, 'gz', gz, 'ofname', 'readout.mod', 'desc', 'DESS stack-of-spirals readout');
T = 4e-3*[0.5:1:rhfrsize];
figure; plot(T,gx,'r'); hold on; plot(T,gy,'g'); plot(T,gz,'b');
xlabel('time (msec)'); ylabel('G/cm');



%% create scan,dess,spiral.tgz
writeloop;

fprintf(1,'\n\tcreate tar archive...');
system('./tarit');
fprintf(1,'\tdone\n');

fprintf(1,'Calculating sequence for display...');
d = readloop('scanloop.txt');
plotseq(800, 820, 'loopArr', d);
%system(sprintf('cp scan.tgz scan_te%dms_tr%.1fs.tgz',round(TE/1e3),TR/1e6)');
fprintf(1,'\n\tdone\n');

fprintf('Next step: Copy scan,dess,spiral.tgz to /usr/g/bin/ and extract, then scan with toppev2\n');


