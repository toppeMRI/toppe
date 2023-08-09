%setup;

addpath('./helper_funcs')
sc = true;    % sc for "scan", set to true when actually running the scan 
              % and you want to send the files to the scanner
 
%%
seq01_Fid;
pulsegeq.seq2ge(get_cur_seqfile);
update_modules(1,2) % module1.mod => tipdown.mod, module2.mod => readout.mod, 
                    % Updates modules.txt accordingly
if sc, toppe.utils.coppe('target','outside'); end

%%
seq02_Fid_multipleFAs;
pulsegeq.seq2ge(get_cur_seqfile);
update_modules(1,2)
if sc, toppe.utils.coppe('target','outside'); end

%%
seq03_SE;
seq03a_SE_withBackGrad;
pulsegeq.seq2ge(get_cur_seqfile);
update_modules(1,2)
if sc, toppe.utils.coppe('target','outside'); end

%%
seq04_SE_withSpolers_nodelay;
pulsegeq.seq2ge(get_cur_seqfile);
update_modules(1,4) 
if sc, toppe.utils.coppe('target','outside'); end

%%
seq05_RadialSE_nodelay;
pulsegeq.seq2ge(get_cur_seqfile);
update_modules(1,7) 
if sc, toppe.utils.coppe('target','outside'); end
toppe.playseq(5,'drawpause',false);


%%  --------------------    data conversion   -----------------------  %%

% constants
ps_sr = .000010;        % pulseq sampling rate
tp_sr = .000004;        % toppe sampling rate


%% radial spin echo sequence

% get most recent p-file
data_file = get_cur_pfile;

% load p-file and reshape
data_pfile = toppe.utils.loadpfile(data_file);
[npt, nch, ~, ~, nevents] = size(data_pfile);
data_pfile = reshape(data_pfile,[npt, nch, nevents]); % should be single sli, echo

% crop based on GE sampling rate of 4 us
n = 256;  % todo! update this to get from .seq file!
pulseq_time_seconds = ps_sr * n;
npts_ge = round(pulseq_time_seconds / tp_sr);
data_pfile = data_pfile(1:npts_ge,:,:);

% interpolate to 10 us
ps_time = 0:ps_sr:pulseq_time_seconds-ps_sr;
tp_time = 0:tp_sr:tp_sr*npts_ge-tp_sr;
data_unsorted = zeros(n,1,n); 
for ii = 1:nevents
    data_unsorted(:,1,ii) = interp1(tp_time,data_pfile(:,:,ii),ps_time);
end

data_file_mat = [data_file(1:end-2),'.mat'];
save(data_file_mat, 'data_unsorted')



%% I guess that's it? There are some other ones in that folder:

seq06_gre_live_demo_step0;
pulsegeq.seq2ge('DEMO_gre_step0.seq');
% rename:
% module1.mod => tipdown.mod
% module3.mod => readout.mod
% Update modules.txt accordingly
toppe.playseq(3);

seq07_gre_live_demo;
pulsegeq.seq2ge('DEMO_gre5.seq'); 
% rename:
% module1.mod => tipdown.mod
% module6.mod => readout.mod
% Update modules.txt accordingly
toppe.playseq(4);


return;

sys = toppe.systemspecs;

B = [1 2 3 4 42 44 48 298 300];

for ib = B
	blk = seq.getBlock(ib);
	figure;
	pulsegeq.plotblock(blk, sys.gamma, sys.raster); title(num2str(ib));
end


