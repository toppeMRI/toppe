% Function to set up the PSD for 3D GRE
% Dinank Gupta, 2024
% University of Michigan

function writeloop_3d_gre(freq_ss,seq,system)

% Inputs:
%   freq_ss     slice select frequency scale
%   seq         sequence information
% system        system information


%% Calculate frequency offsets for each slice

% calculate bandwidth of slice select pulse in Hz
sli_bw = seq.tbw/(1e-3*seq.pulsedur);
%% Adding spacing between adjacent slices:
slispacing_bw90 = round(system.gamma*seq.gamp90);

%^ find out what needs to change here to get the correct slice plan
% find offset of 1st slice
freq_ss_start = (freq_ss - ( (sli_bw+slispacing_bw90) * (seq.nsli-1)));

% create vector with all slice offsets, rounding to nearest Hz
freq_offset_all = -freq_ss_start/2 + round(freq_ss_start:(sli_bw+slispacing_bw90):freq_ss);

if numel(freq_offset_all) ~= seq.nsli
    error('Error calculating slice offsets, check freq_ss parameter.')
end

%% Loop through all timepoints, slices

ver = 6;
dab='off';
seq.n_disdaq=100;
adc_flag=0;

toppe.write2loop('setup',system,'version',ver);
for itp = 0:seq.ntp
    n=0;

    for ishot = 1:seq.n_shots_per_vol
        for nd = 1:seq.n_disdaq % n_disdaq will be set to 1 during adc
            kz = seq.z_scale_shot(ishot,max(1,itp));
            rotmatx = seq.rotmat_shot(:,:,ishot,max(1,itp));
            n=n+1;

            %tipdown
            toppe.write2loop('tipdown.mod',system,'RFoffset',freq_offset_all, ...
                'echo',1,'slice',1,'view',1, 'RFspoil',1,...
                'core',1,'version',ver);
            % Delay module to get the correct TE timing
            toppe.write2loop('delay',system,'textra',seq.delay_TE,'core',1)

            if(seq.adddwgrad)
                % dw grad
                toppe.write2loop('arfigrad.mod',system,...
                    'echo',1,'slice',1,'view',1,...
                    'version',ver,'core',1,...
                    'Gamplitude',[1;1;1]*seq.grad_dir(max(1,itp)));
            end

            % spiral readout
            toppe.write2loop('readout.mod',system,...
                'echo',1,'slice',n,'view',max(1,itp), 'RFspoil',1,...
                'rotmat',rotmatx,'version',ver,...
                'Gamplitude',[1;1;adc_flag*kz],...
                'dabmode',dab, 'core',2);


            toppe.write2loop('crusher.mod',system,...
                'echo',1,'slice',1,'view',1,...
                'version',ver, 'core',3);

            % Delay module to get the correct TR
            toppe.write2loop('delay',system,'textra',seq.delay_TR,'core',3)
        end
        if(adc_flag==0)% After the disdaq loop ends once, set adc to 1.
            adc_flag=1;
            seq.n_disdaq=1;
            dab='on';
            break
        end

    end

end



toppe.write2loop('finish',system);


% Create 'sequence stamp' file for TOPPE.
% This file is listed in line 6 of toppe0.entry
toppe.preflightcheck(['toppe',num2str(seq.cv_num),'.entry'], 'seqstamp.txt', system);
copyfile(['toppe',num2str(seq.cv_num),'.entry'],'toppeN.entry')

end

