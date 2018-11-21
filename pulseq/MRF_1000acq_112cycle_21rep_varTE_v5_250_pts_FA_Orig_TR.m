%%Test making sequence with fixed FA and TR
clear
clc

%% direct to pulseq folder
Pulseq = 'C:\Users\qiane\Desktop\Columbia\pulseq-master-1.1\matlab';
addpath(genpath('.'));
addpath(genpath(Pulseq));
os = 'pc';

%% Set system limits
gamma = 42576000; % in Hz %Determined from Pulseq - do not change,%Hz/T
Gmax = 32; %mT/m
SRmax = 130;%T/m/s
system = mr.opts('MaxGrad',Gmax,'GradUnit','mT/m',...
 'MaxSlew',SRmax,'SlewUnit','T/m/s', 'gradRasterTime', 10e-6, 'rfRingdownTime',10e-6);
seq=mr.Sequence(system);
Cycle_len = 100;
total_acq = 1000;
RepeatTimes = ceil(total_acq/Cycle_len);

%% Load data
fov=225e-3;
Nx=256;% Ny=256;
sliceThickness=5e-3;
load('method_orig.mat')
TR_all = (method.VariableTR+6)*1e-3;%converted to s
% TE = 2e-3;%minimum - changed from 5, when fixed is 2
dx = fov/Nx;
dy = dx;
SRmax = SRmax.*0.8;

%% setting up trajectory
Nshots = 48; %# of spirals
alpha = 1;%constant spiral density
[ktraj,dcf,t,ind,out,G]=design_spiral_pulseq(fov*1000,Nx,Nshots,10,...
    '',Gmax,SRmax,'1H',true,true,true);%fov is in mm
G = G.*gamma; %converting T/m to Hz/m

%% Rotating G, the G from design_spiral_pulseq function has rewinders but no rotation
% ktraj has rotations but no rewinders
G = complex(G(1,:),G(2,:));
G = [0,G];
Gn = zeros(size(G,2),Nshots);
for ks = 1:Nshots
    Gn(:,ks) = (G.*exp(1i*(ks-1)/Nshots*2*pi))';
end

%% Break ktraj into Nshots
kshot = zeros(length(ktraj)/Nshots,Nshots);
OneSpiralPoints = length(ktraj)/Nshots; %# of points in one spiral
for n1 = 1:Nshots
    kshot(:,n1) = ktraj((n1-1)*OneSpiralPoints+1:n1*OneSpiralPoints).*1e3;% per m
    plot (kshot(:,n1));hold on;grid on;xlabel('kx (m-1)'); ylabel('ky(m-1)');
end

%% Get gradient based on k data
DT=system.gradRasterTime;
ktrajs = zeros(size(ktraj,1), size(ktraj,2), 2);
ktrajs(:,:,1) = real(ktraj);
ktrajs(:,:,2) = imag(ktraj);

%% Slice selection, also dowmsample the FA data so that its peroid is total_acq
adc = mr.makeAdc(length(Gn),'Dwell',system.gradRasterTime);
FA_shots = GenerateFA;

for itr = 1:100
    temp1 = (itr-1)*10+1;
    [rf, gz] = mr.makeSincPulse(deg2rad(FA_shots(temp1)),system,'Duration', 2e-3,...
        'SliceThickness',sliceThickness,'apodization',0.5,'timeBwProduct',4); %for GE RF pulse duration needs to be constant
    rf_all(itr) = rf; 
    dur_rf_all(itr) = mr.calcDuration(rf_all(itr));
    gz_all(itr) = gz;
end

%% Need to determine multi-slice frequency offsets
Nslices=1;
deltaz=Nslices*sliceThickness;
z=-(deltaz/2):sliceThickness:(deltaz/2);
gzReph = mr.makeTrapezoid('z',system,'Area',-gz.area/2,'Duration',1e-3);
gx = mr.makeArbitraryGrad('x', squeeze(real(Gn(:,1))),system);
[rf2] = mr.makeBlockPulse(pi,system,'Duration',500e-6);
rf2.deadTime = 100e-6;rf2.ringdownTime=30e-6;%changing properties only for the rf block
gzSpoil = mr.makeTrapezoid('z',system,'Area',gz.area*2,'Duration',3*8e-4);

%% Define sequence block
numfullim = round(length(TR_all)/Nshots);
delay18 = 18e-3-mr.calcDuration(gzSpoil)-mr.calcDuration(rf2)/2;
delaySpoil = mr.makeDelay(delay18);

seq.addBlock(gzSpoil);
seq.addBlock(rf2);
seq.addBlock(gzSpoil);
seq.addBlock(delaySpoil);
TE = mr.calcDuration(gzReph) +max(dur_rf_all) + 0.1e-3; %minTE
time =0;
for ds = 1:numfullim
    for ns=1:Nshots
        
        n = (ds-1)*Nshots+ns;
        if n <= 1000
            
            n2 = mod(n,Cycle_len); %to avoid rf pulse shape issues
            if(n2==0)
                n2 =Cycle_len;
            end
            TR = TR_all(n);
            %% Calculate timing
            %         delayTE= TE - mr.calcDuration(gzReph) - (mr.calcDuration(rf_all(n2))/2);
            delayTE = 0.5e-3;
            TE = mr.calcDuration(gzReph) + (mr.calcDuration(rf_all(n2))/2)+ delayTE;
            TE_all(n2) = TE;
            % %         delayTE = delayTE - mod(delayTE, 1e-5);%check back later5
            if(TR < (TE + mr.calcDuration(gx)))
                TR = TE + mr.calcDuration(gx) +0.1e-3;
                TR_all(n) = TR;
            end
            
            delayTR= TR - TE - mr.calcDuration(gx);
            delayTR = delayTR - mod(delayTR, 1e-5);
            
            delay1 = mr.makeDelay(delayTE);
            delay2 = mr.makeDelay(delayTR);
            
            seq.addBlock(rf_all(n2),gz_all(n2));
            seq.addBlock(gzReph);
            gx = mr.makeArbitraryGrad('x', squeeze(real(Gn(:,ns))),system);
            gy = mr.makeArbitraryGrad('y', squeeze(imag(Gn(:,ns))),system);
            seq.addBlock(delay1);
            seq.addBlock(gx,gy,adc);
            seq.addBlock(delay2);
            time = time + TR;
        end
    end
end
% disp(time/60);
% seq.plot('TimeRange',[0 TR]);
fname = ['Spiral_2D_vari_1000_single_FAScale=1_4GE_', num2str(Nshots),'_',num2str(TE)];
% fname = ['7'];
seq.write([fname,'.seq']);
% scan_time = sum(TR_all).*Nslices/60;
% disp(scan_time);
