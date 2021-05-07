% very basic and crude non-Cartesian recon using griddata()
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
data_file_path=[path D(I(end-0)).name] % use end-1 to reconstruct the second-last data set, etc.

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

[ktraj_adc_nom,t_adc] = seq.calculateKspacePP('trajectory_delay',0e-6); 

% detect slice dimension
max_abs_ktraj_adc=max(abs(ktraj_adc_nom'));
[~, slcDim]=min(max_abs_ktraj_adc);
encDim=find([1 2 3]~=slcDim);

ktraj_adc_nom = reshape(ktraj_adc_nom, [3, adc_len, size(ktraj_adc_nom,2)/adc_len]);

prg_angle=squeeze(atan2(ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:),ktraj_adc_nom(encDim(2),2,:)-ktraj_adc_nom(encDim(2),1,:)));
nproj=length(prg_angle);

i_pure2=find(abs(ktraj_adc_nom(encDim(2),1,:))<eps);
i_pure1=find(abs(ktraj_adc_nom(encDim(1),1,:))<eps);
assert(length(i_pure2)==length(i_pure1)); % the code below assumes this

delta180=nproj/length(i_pure2);

%%

data_fft1=ifftshift(ifft(ifftshift(rawdata,1)),1);

ip=1:nproj;
cmplx_diff=data_fft1(2:end,ip,:).*conj(data_fft1(end:-1:2,mod(ip-1+delta180,nproj)+1,:));
cmplx_diff_no_channels=sum(cmplx_diff,3);

figure; imagesc(angle(cmplx_diff_no_channels));

%% pick the pure X and pure Y differences, plot them and extimate the slope

thresh=0.05;
cmplx_diff_pure_axes=cmplx_diff_no_channels(:,[i_pure2 i_pure1]);
mpa=max(abs(cmplx_diff_pure_axes));
for i=1:size(cmplx_diff_pure_axes,2)
    cmplx_diff_pure_axes(abs(cmplx_diff_pure_axes(:,i))<thresh*mpa(i),i)=0;
end

% plot
figure; plot(angle(cmplx_diff_pure_axes(:,[1 3])));

% just get the mean slope of the phase
msop=angle(sum(cmplx_diff_pure_axes(2:end,:).*conj(cmplx_diff_pure_axes(1:end-1,:))));
delays_imageSpace=-msop/2/pi/2*(t_adc(2)-t_adc(1))*adc_len; 

%% a more accurate function
delays=[0 0 0];
delays(encDim)=calc_delays(rawdata, ktraj_adc_nom(encDim,:),i_pure2,i_pure1,t_adc(2)-t_adc(1));

% improve iteratively & verify
n_it=4;
d=[0 0 0];
for n=1:n_it
    d=sum(delays,1);
    ktraj_adc = seq.calculateKspacePP('trajectory_delay',d);
    d(encDim)=calc_delays(rawdata, ktraj_adc(encDim,:),i_pure2,i_pure1,t_adc(2)-t_adc(1));
    delays = [delays; d];
end
d=sum(delays,1);
fprintf('found delays: [%g %g %g] us\n', round(d*1e9)/1e3);

%%

function delays=calc_delays(rawdata, ktraj_adc,i_pureX,i_pureY,dt)
    adc_len=size(rawdata,1);
    n_chan=size(rawdata,3);
    n_sel=length(i_pureX)+length(i_pureY);
    data_selected=rawdata(:,[i_pureX i_pureY],:);
    data_sel_res=zeros(adc_len,n_sel,n_chan);
    dk=zeros(1,n_sel);abs(diff(ktraj_adc(1,(i_pureX(1)-1)*adc_len+(1:2))));
    for i=1:n_sel
        if i<=length(i_pureX)
            kr_selected=ktraj_adc(1,(i_pureX(i)-1)*adc_len+(1:adc_len));
        else
            kr_selected=ktraj_adc(2,(i_pureY(i-length(i_pureX))-1)*adc_len+(1:adc_len));
        end
        dk(i)=diff(kr_selected(1:2));
        for c=1:n_chan
            data_sel_res(:,i,c)=interp1(kr_selected,data_selected(:,i,c),(-adc_len/2:(adc_len/2-1))*abs(dk(i)),'pchip',0);
        end
    end

    data_sel_fft1=ifftshift(ifft(ifftshift(data_sel_res,1)),1);

    thresh=0.05;
    cmplx_diff_sel=data_sel_fft1(:,2:2:end,:).*conj(data_sel_fft1(:,1:2:end,:));
    cmplx_diff_sel_no_channels=squeeze(sum(cmplx_diff_sel,3));
    mpa=max(abs(cmplx_diff_sel_no_channels));
    for i=1:size(cmplx_diff_sel_no_channels,2)
        cmplx_diff_sel_no_channels(abs(cmplx_diff_sel_no_channels(:,i))<thresh*mpa(i),i)=0;
    end
    
    figure; plot(angle(cmplx_diff_sel_no_channels));

    % just get the mean slope of the phase
    msop1=angle(sum(cmplx_diff_sel_no_channels(2:end,:).*conj(cmplx_diff_sel_no_channels(1:end-1,:))));
    delays=msop1/2/pi/2*dt*adc_len.*sign(dk(1:2:end)); % we need thins "sign" to account for the direction (positive or negative) of the first projection in each pair
end
