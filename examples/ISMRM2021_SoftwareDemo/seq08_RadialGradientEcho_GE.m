% setup
seq08_RadialGradientEcho;             % create 'gre_rad.seq'
pulsegeq.seq2ge('gre_rad.seq');       % create TOPPE files (modules.txt, scanloop.txt, *.mod files)
nModulesPerTR = 4;                    % by inspection
toppe.playseq(nModulesPerTR);         % preview sequence (loop)

