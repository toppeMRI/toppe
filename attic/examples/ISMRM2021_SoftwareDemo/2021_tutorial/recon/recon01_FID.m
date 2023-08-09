% very basic FID data handling example
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


%% we just want time stamps: t_adc, t_excitation
[~,t_adc,~,~,t_excitation]=seq.calculateKspacePP();

t_excitation_adc=t_excitation((end-readouts+1):end); % remove dummy/prep scans
t_adc=reshape(t_adc,[adc_len,readouts]);
t_e=t_adc-t_excitation_adc(ones(1,adc_len),:);

% plot raw data
figure; plot(t_e, abs(rawdata(:,:))); title('raw signal(s)'); xlabel('time /s');
hold on; plot(t_e(1),0); % trick to force Y axis scaling

if (readouts>1)
    figure; plot(abs(rawdata(1,:))); title('time evolution (first FID point)'); 
    dataavg=squeeze(mean(rawdata,2));
    figure; plot(t_e, abs(dataavg)); title('averaged FIDs');
end

% calculate frequency axis
f=zeros(size(t_e));
adc_dur=(t_adc(end,1)-t_adc(1,1))/(adc_len-1)*adc_len;
f(:)=1/adc_dur;
f=cumsum(f,1)-0.5/adc_dur*adc_len;
% plot spectra
figure; plot(f, abs(ifftshift(ifft(ifftshift(rawdata(:,:),1)),1))); title('spectr(um/a)'); xlabel('frequency /Hz');
