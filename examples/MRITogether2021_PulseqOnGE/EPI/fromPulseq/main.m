% Convert Pulseq file created by ../epi.m back to TOPPE.

addpath ~/github/pulseq/matlab/       % +mr package
addpath ~/github/toppeMRI/PulseGEq/   % +pulsegeq package

% Create TOPPE files
pulsegeq.seq2ge('epi.seq', 'toppeVersion', 'v4', 'system', seq.sys);

% Display scan loop
figure;
iStart = 1;  % row in scanloop.txt
iStop = 10;
toppe.plotseq(iStart, iStop, seq.sys);
