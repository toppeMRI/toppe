addpath ~/github/pulseq/              % +mr package
addpath ~/github/toppeMRI/PulseGEq/   % +pulsegeq package
addpath ~/github/toppeMRI/toppe/      % +toppe package
seq08_RadialGradientEcho;             % create 'gre_rad.seq'
pulsegeq.seq2ge('gre_rad.seq');       % create TOPPE files (modules.txt, scanloop.txt, *.mod files)
toppe.playseq(4);                     % preview sequence (loop)

