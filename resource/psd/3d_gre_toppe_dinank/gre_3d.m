% 3D GRE using stack of spirals
% Coded in  toppe v6
% Dinank Gupta, 2024
% University of Michigan


%% Path


addpath(genpath('/home/dinankg/tools/toppev6')) %get tv6
clear
clc


import toppe.*
import toppe.utils.*

%%
ACTUAL_SCAN = 1
seq.scanner = 'inside';

seq.FOV = 24;           % Starting Field of view, x-y, cm
seq.nshot = 4;         % For multishot spiral
seq.npix = 128;          % Matrix size, x-y
seq.und_samp_fact = 1 % Undersampling factor. Should be divisble by nshot for now.

seq.cv_num = str2double(sprintf('%d', 72, seq.FOV)); %enrty file number.

%% Setup sequence parameters using getparams function

% KEY parameters set here

seq.tmp_OSSI_sli_offset =0; %4.89;%2.55;  %% volume offset from isocenter - in cm!!

seq.zFOV = 4.5;          % z FOV, cm
seq.nsli = 1;            % Number of slices
seq.zsteps=15; % kz encodes

seq.ntp = 5;
seq.spiraldir = 2;     % 1 is spiral out, 2 is spiral in
sp_sl_sc = .7;  % vds spiral slew scaling 0.7
% crash
seq.dt = 4e-6;          % TOPPE sampling rate in seconds


% most parameters set here (one's that aren't changed as often)
seq.slthick = seq.zFOV;     % Slice thickness in cm %for 3D excitation
% Spatial settings


seq.TR = 39;
seq.TE = 21;

%% flags:
seq.adddwgrad=0; % To have dw grad or not.
% Using Ernst angle for T1 of 1500ms
seq.fa=acosd(exp(-seq.TR/1500));

%% Option for real scan or dev
%  "ACTUAL_SCAN" runs toppe.coppe and saves .mod files

if ACTUAL_SCAN
    exp_number = 0;
    while true
        exp_number = exp_number + 1;
        exp_name = ['experiments/',datestr(now,'yyyy-mm-dd_'),seq.scanner,'_exp',num2str(exp_number)];
        if ~isfolder(exp_name)
            mkdir(exp_name)
            break
        end
    end
end



%% setup toppe system using defaults

% system = toppe.systemspecs();
system = toppe.systemspecs('maxSlew', 20,'maxRF',0.20);
ncycles = 0; % balanced sl select.
if(strcmp(seq.scanner, 'inside'))
    system.maxGrad = 10;
else
    system.maxGrad = 5;
end

%% make tipdown.mod
seq.pulsedur = 4;      % RF pulse duration in ms
seq.tbw = 6;           % Time bandwidth for SLR pulse
seq.rf_sl_sc = 0.76;    % rf slew scale for calling makeslr function

if seq.fa < 20, type = 'st'; else, type = 'ex'; end
freq_ss=0;
[rf, gex, freq_ss, fnamestem] = toppe.utils.rf.makeslr(seq.fa, seq.slthick, seq.tbw, ...
    seq.pulsedur, eps, system, ...
    'maxSlewScale',seq.rf_sl_sc,'type',type,'ofname', 'tipdown.mod');

rf = cat(1,rf,zeros(25,1));
gex = cat(1,gex,zeros(25,1));

seq.gamp90 = max(gex);

toppe.writemod(system, 'rf', rf, 'gz', gex,'nChop',[25,0], 'ofname', 'tipdown.mod');


%% OPTIONAL: Adding bipolar ARFI grads to this
dwamp = 4; %G/cm
dwlen = 4 * 1e-3; %sec
dw_dir=[0 1 0];
dw_gradient = toppe.utils.trapwave2(dwamp*dwlen, dwamp, 7, seq.dt*1e3);
dw_gx=dw_dir(1)*dw_gradient';
dw_gy=dw_dir(2)*dw_gradient';
dw_gz=dw_dir(3)*dw_gradient';

% Delay between the 2 gradients to allow ultrasound to rise to max
seq.dw_grad_delay = 1; %ms
grad_delay_samp = seq.dw_grad_delay /seq.dt/1e3;
seq.num_dw_grad = 2; % For OGSE type SE.
% Flip and append
dw_gxall = dw_gx; dw_gyall = dw_gy; dw_gzall = dw_gz;
for g =1:seq.num_dw_grad - 1
    dw_gxall = [dw_gxall;zeros(grad_delay_samp,1);dw_gx(2:end)*(-1)^g];
    dw_gyall = [dw_gyall;zeros(grad_delay_samp,1);dw_gy(2:end)*(-1)^g];
    dw_gzall = [dw_gzall;zeros(grad_delay_samp,1);dw_gz(2:end)*(-1)^g];
end
% Time to send the trigger (us):
quick_trig=0;
if(seq.num_dw_grad==1 && quick_trig==1)
    t_trigout = 200;
elseif(seq.num_dw_grad==1)
    t_trigout = 1000;
else
    t_trigout = length(dw_gradient)* seq.dt*1e6+100;
end
seq.grad_dir = ones(seq.ntp,1);
toppe.writemod(system,'gx',dw_gxall,'gy',dw_gyall,'gz',dw_gzall,'ofname','arfigrad.mod');
%% Making list of triggers
seq.triglist = ones(1,seq.ntp);
% Turning off
seq.triglist(1,[1,4]) = 0;


%% Generate balanced VD spiral, WITH CRUSHER for sens maps

seq.oversamp = 600;

[gx,gy,t,npts] = makevdspiral2(seq.FOV,seq.npix,seq.nshot,seq.oversamp,5,sp_sl_sc*system.maxSlew*10);
paramsint16(2) = npts;% p
gz=0*gx;
gx=gx(:,1);gy=gy(:,1);gz=gz(:,1);
if(seq.spiraldir==2)
    fprintf('Not tested with spiral in')
    gx=flip(gx);gy=flip(gy);gz=flip(gz);
end
% List of rotations for spirals:
shotrot = linspace(0,2*pi-(2*pi)/seq.nshot*seq.und_samp_fact,seq.nshot/seq.und_samp_fact);
seq.rot_list =cat(2,repmat([0 0 1],seq.nshot/seq.und_samp_fact,1),shotrot.');
seq.rotmat = axang2rotm(seq.rot_list);

%% Make steps in z
gambar = 4257;                           % gamma/2pi in Hz/Gauss
dkz = 1/seq.zFOV; %1/cm
res = seq.zFOV/seq.zsteps;%cm
seq.kzmax = 1/res;%1/cm
Gz_area = seq.kzmax/gambar; %G sec/cm grad area to reach kz max

% Grads will go from -Gz_area/2 to Gz_area/2
zp = (-(seq.zsteps-1)/2:(seq.zsteps-1)/2)/((seq.zsteps-1)/2)/2; %zpe scaling for each TR
%List of gradient area:
gz_trap_area_list = -1*Gz_area; %G sec/cm
%First trap which will be the longest one.
z_enc_trap = toppe.utils.trapwave2(gz_trap_area_list(1),6, system.maxSlew-5,system.raster*1e-3).';
seq.zp_scale = zp; %Scaling of the kz traps
% Appending z enc to spirals.
% z_enc_trap = repmat(z_enc_trap,1,seq.nshot);
gz=cat(1,z_enc_trap,gz,-1*z_enc_trap);
gx=cat(1,0*z_enc_trap,gx,0*z_enc_trap);
gy=cat(1,0*z_enc_trap,gy,0*z_enc_trap);
paramsint16(14) = length(z_enc_trap);
%NOTE than total number of samples to skip will be paramsint16(14)+paramsint16(11)
% And the spiral in samples are paramsint16(2) for recon.
toppe.writemod(system,'gx',gx,'gy',gy,'gz',gz,'ofname','readout.mod','hdrints',paramsint16);

%% Make a table of rotation matrices and kz-encode scaling.
GR = 1.618;
seq.doGA=0; % Flag to turn GA on/off
GA_tp_n = randi(1e2,seq.ntp,1); % GA differences across time.
seq.rot_angle_shot=[];seq.z_scale_shot=[];% shots x tp
for nt=1:seq.ntp
    seq.n_shots_per_vol=0;
    for nspiral = 1:(seq.nshot/seq.und_samp_fact)
        n=0; %Resetting the golden angle rotations for the same kz platters,
        %just adding rotations
        for nkz = 1:seq.zsteps
            seq.n_shots_per_vol=seq.n_shots_per_vol+1;
            n=n+1;
            seq.rot_angle_shot(seq.n_shots_per_vol,nt) = mod(seq.doGA*(360-(360/GR))*(n+GA_tp_n(nt))+...
                rad2deg(shotrot(nspiral)),360);
            % seq.rot_angle_shot(seq.n_shots_per_vol,nt) = 360/nspiral;
            seq.z_scale_shot(seq.n_shots_per_vol,nt) = seq.zp_scale(nkz);
        end
    end
end
seq.rot_angle_shot = deg2rad(seq.rot_angle_shot);
seq.rot_list_shot =cat(2,repmat([0 0 1],seq.n_shots_per_vol*seq.ntp,1),seq.rot_angle_shot(:));
seq.rotmat_shot = axang2rotm(seq.rot_list_shot);
seq.rot_list_shot = reshape(seq.rot_list_shot,[],seq.ntp,4);%shots x tp x 4
seq.rotmat_shot = reshape(seq.rotmat_shot,3,3,[],seq.ntp);% 3 x3 x shots x tp
%% Crusher at the end of spiral.
gspoil = toppe.utils.makecrusher(4, seq.slthick/seq.zsteps, system, 0, 8,6);
toppe.writemod(system,'gx',gspoil,'gy',gspoil,'gz',gspoil,'ofname','crusher.mod','hdrints',paramsint16);
%% Getting min TE and appending time to get the right TE
n_rf_left = length(rf) - (seq.pulsedur/2)/seq.dt/1e3;

minTE = (seq.adddwgrad*length(dw_gxall) + n_rf_left + ...
    length(gx)*(seq.spiraldir == 2))*seq.dt*1e3;
if(minTE>seq.TE)
    error(['TE is too small, set it to atleat ',num2str(minTE),' ms'])
end
seq.delay_TE  = max(0,seq.TE - minTE);

%% Making modules and cores.txt
% Entries are tab-separated.
% Calculating the min durations for each mod:
tipdown_dur = 100+length(rf)*seq.dt*1e6;
readout_dur = length(gx)*seq.dt*1e6;
crusher_dur = length(gspoil)*seq.dt*1e6;
arfi_dur = length(dw_gxall)*seq.dt*1e6;
% Entries are tab-separated.
modFileText = ['' ...
    'Total number of unique cores\n' ...
    '4\n' ...
    'fname	duration(us)	hasRF?	hasDAQ? hastrigout?\n' ...
    'tipdown.mod	',num2str(tipdown_dur),'	1	0   -1\n'...
    'readout.mod	',num2str(readout_dur),'	0	1   -1\n' ...
    'crusher.mod	',num2str(crusher_dur),'	0	0   -1\n'...
    'arfigrad.mod ',num2str(arfi_dur),'   0   0   -1'];
fid = fopen('modules.txt', 'wt');
fprintf(fid, modFileText);
fclose(fid);

% To-do : add wait to acquire FMAPS.
if(seq.adddwgrad)
    coreFileText{1} = [1,4,0];%tip-grad-<delay>
    coreFileText{2} = [2];%kz-spiral-kz
    coreFileText{3} = [3,0];%crush-<delay>
else
    coreFileText{1} = [1,0];%tip-<delay>
    coreFileText{2} = [2];%-kz-spiral-kz
    coreFileText{3} = [3,0];%crush-<delay>
end

writecoresfile(coreFileText);

%% Setting TR
time_current = (seq.adddwgrad*arfi_dur + tipdown_dur + readout_dur + crusher_dur)*1e-3+seq.delay_TE;
time_add = seq.TR - time_current;
if(time_add<0)
    error('Increase TR')
end

seq.delay_TR= time_add;

writeloop_3d_gre(freq_ss, seq, system);

%% Copy to scanner and save TOPPE files

if ACTUAL_SCAN
    [status] = copyfile('readout.mod', ['./',exp_name,'/readout.mod']);
    [status] = copyfile('tipdown.mod', ['./',exp_name,'/tipdown.mod']);
    [status] = copyfile('arfigrad.mod', ['./',exp_name,'/arfigrad.mod']);
    [status] = copyfile('crusher.mod', ['./',exp_name,'/crusher.mod']);

    %     [status] = copyfile('tipdown_fatsat.mod', ['./',exp_name,'/tipdown_fatsat.mod']);
    [status] = copyfile('scanloop.txt', ['./',exp_name,'/scanloop.txt']);
    [status] = copyfile('seqstamp.txt', ['./',exp_name,'/seqstamp.txt']);

    [status] = copyfile('toppe-scanfiles.tgz', ['./',exp_name,'/toppe-scanfiles.tgz']);
    [status] = copyfile('modules.txt', ['./',exp_name,'/modules.txt']);
    save('seq.mat','seq')
    [status] = copyfile('seq.mat', ['./',exp_name,'/seq.mat']);
    toppe.utils.coppe('target',seq.scanner,'cv',seq.cv_num,'use_pw',true,'version',6);
end

%% Sequence time:
time_msec = seq.TR*seq.zsteps*seq.nshot
%% Play sequence out virtually
figure;

for j=1:11
    toppe.plotseq(1+5*(j-1),5*j,system,'doDisplay','true','gmax',system.maxGrad);
    drawnow
end
if ACTUAL_SCAN
    [status] = copyfile('psd.png', ['./',exp_name,'/psd.png']);
end

