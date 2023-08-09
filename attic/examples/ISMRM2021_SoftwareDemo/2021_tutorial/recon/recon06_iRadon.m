% very basic inverse Radon radial reconstruction
%
% it loads Matlab .mat files with the rawdata in the format 
%     adclen x channels x readouts
% it also seeks an accompanzing .seq file with the same name to interpret
%     the data

%% Load the latest file from a dir
path='../data_ocra/'; % directory to be scanned for data files
%path='../IceNIH_RawSend/demo-set/'; % alternative
pattern='*.mat';
D=dir([path pattern]);
[~,I]=sort([D(:).datenum]);
data_file_path=[path D(I(end-2)).name] % use end-1 to reconstruct the second-last data set, etc.

data_unsorted = load(data_file_path);

if isstruct(data_unsorted)
    fn=fieldnames(data_unsorted);
    assert(length(fn)==1); % we only expect a single variable
    data_unsorted=data_unsorted.(fn{1});
end

%% Load sequence from file 
seq = mr.Sequence();              % Create a new sequence object
seq_file_path = [data_file_path(1:end-3) 'seq'];
%seq.read(seq_file_path); % no because it conflicts with the FID multiple FAs
seq.read(seq_file_path,'detectRFuse'); % this is a better option for the SE sequences

%% raw data preparation

if ndims(data_unsorted)<3 && size(data_unsorted,1)==1
    % OCRA data may need some fixing
    [~, ~, eventCount]=seq.duration();
    readouts=eventCount(6);
    data_unsorted=reshape(data_unsorted,[length(data_unsorted)/readouts,1,readouts]);
end

[adc_len,channels,readouts]=size(data_unsorted);
rawdata = permute(data_unsorted, [1,3,2]); % channels last

%% Analyze the nominal trajectory

[ktraj_adc_nom,t_adc] = seq.calculateKspacePP('trajectory_delay',0); 

% detect slice dimension
max_abs_ktraj_adc=max(abs(ktraj_adc_nom'));
[~, slcDim]=min(max_abs_ktraj_adc);
encDim=find([1 2 3]~=slcDim);

ktraj_adc_nom = reshape(ktraj_adc_nom, [3, adc_len, size(ktraj_adc_nom,2)/adc_len]);

prg_angle=unwrap(squeeze(atan2(ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:),ktraj_adc_nom(encDim(1),2,:)-ktraj_adc_nom(encDim(1),1,:))));
nproj=length(prg_angle);

%% from k-space to projections (1D FFTs)

data_fft1=ifftshift(ifft(ifftshift(rawdata,1)),1);

figure; imagesc(abs(squeeze(data_fft1))); title('sinogramm view')

%% crop the data to remove oversampling and adapt to the iradon
target_matrix_size=200;
shift=-3; % 12 was found experimentally for the first Benjamin's data set
cropLeft=(adc_len-target_matrix_size)/2+shift;
cropRight=(adc_len-target_matrix_size)/2-shift;
data_fft1c=data_fft1(2+cropLeft:end-cropRight,:,:);

% visualize the matching of positive and negative directions 
p1=1;
[~,p2]=min(abs(mod(prg_angle-(prg_angle(p1)-pi)+pi,2*pi)-pi)); % look for the projection closest to the opposite to p1 
figure;plot(abs(data_fft1c(1:end,p1,1)));hold on;plot(abs(data_fft1c(end:-1:1,p2,1))); title('comparing opposite projections');

%% the actuall iRadon transform
theta=180-prg_angle/pi*180;
for c=1:channels
    % the classical (absolute value) transform
    irad_a=iradon(abs(data_fft1c(:,:,c)),theta,'linear','Hann');
    %irad_a=iradon(abs(data_fft1c(:,:,c)),theta);
    %irad_a=iradon(abs(data_fft1c(:,:,c)),theta,'linear','Shepp-Logan');
    if (c==1)
        irad_abs=zeros([size(irad_a) channels]);
        irad_cmpx=zeros([size(irad_a) channels]);
    end
    irad_abs(:,:,c)=irad_a;
    % MR-specific complex-valued transform
    irad_r=iradon(real(data_fft1c(:,:,c)),theta,'linear','Hann');
    irad_i=iradon(imag(data_fft1c(:,:,c)),theta,'linear','Hann');
    irad_cmpx(:,:,c)=irad_r + 1i*irad_i;
end

figure;imab(abs(irad_abs));colormap('gray'); title('abs iRadon recon')
figure;imab(abs(irad_cmpx));colormap('gray'); title('complex iRadon recon')
%axis('equal');

%% Sum of squares combination
if channels>1
    sos=abs(sum(irad_cmpx.^2,ndims(irad_cmpx)).^(1/2));
    sos=sos./max(sos(:));
    figure;imab(sos);colormap('gray'); title('SOS complex iRadon recon')
    %imwrite(sos, ['img_combined.png'])
end

return 

%% try to predict projections based on the iRadon recon and correct the phase of the actual measured data with respect to the syntheric projection data 

rad_abs = radon(irad_abs,theta);
rad_abs=rad_abs(2:end-1,:); % I assume the projecttions have grown because of the roinding errors

rad_abs_fft=fftshift(fft(fftshift(rad_abs,1)),1);
data_abs_fft=fftshift(fft(fftshift(abs(data_fft1c),1)),1);

%% try to detect frequency shifts on the simulated k-space
figure;imagesc(angle(data_abs_fft.*conj(rad_abs_fft)));

%figure;plot(angle(data_abs_fft.*conj(rad_abs_fft)));
%figure;plot(abs(data_abs_fft.*conj(rad_abs_fft)));

data_rad_diff=data_abs_fft.*conj(rad_abs_fft); 
mx=max(abs(data_rad_diff(:)));
data_rad_diff(abs(data_rad_diff(:))<0.02*mx)=0;

shift=angle(sum(data_rad_diff(2:end,:).*conj(data_rad_diff(1:end-1,:))))/2/pi*size(rad_abs_fft,1);
figure;plot(shift); title('detected shifts in pixels');

%% try to look at the phase slopes and the mean phase in the image domain 

diff_ra=data_fft1c.*rad_abs;
mx=max(abs(diff_ra(:)));
diff_ra(abs(diff_ra(:))<0.1*mx)=0;

slope=angle(sum(diff_ra(2:end,:).*conj(diff_ra(1:end-1,:))));
rr=-(size(diff_ra,1)-1)/2:(size(diff_ra,1)-1)/2;

diff_ra_noslope=diff_ra.*exp(-1i*rr'*slope);

avg_ph=angle(sum(diff_ra_noslope));

%figure; imagesc(angle(diff_ra));
%figure; plot(angle(diff_ra_noslope));
%figure; plot(slope);

%% phase-correct the data

data_fft1c_pc=data_fft1c.*exp(-1i*(rr'*slope+avg_ph));

for c=1:channels
    % MR-specific complex-valued transform
    irad_r=iradon(real(data_fft1c_pc(:,:,c)),theta,'linear','Hann');
    irad_i=iradon(imag(data_fft1c_pc(:,:,c)),theta,'linear','Hann');
    if (c==1)
        irad_cmpx_pc=zeros([size(irad_r) channels]);
    end
    irad_cmpx_pc(:,:,c)=irad_r + 1i*irad_i;
end

figure;imab(abs(irad_cmpx_pc));colormap('gray'); title('complex phase-corrected iRadon recon')
