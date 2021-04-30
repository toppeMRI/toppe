%setup;

seq01_Fid;
pulsegeq.seq2ge('fid.seq');
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly

seq02_Fid_multipleFAs;
pulsegeq.seq2ge('fid.seq');
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly

seq03_SE;
seq03a_SE_withBackGrad;
pulsegeq.seq2ge('se.seq');
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly

seq04_SE_withSpolers;
pulsegeq.seq2ge('se.seq');
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly

seq05_RadialSE;
pulsegeq.seq2ge('se_radial.seq');
% rename:
% module1.mod => tipdown.mod
% module6.mod => readout.mod
% Update modules.txt accordingly
toppe.playseq(4,'drawpause',false);

seq06_gre_live_demo_step0;
pulsegeq.seq2ge('DEMO_gre_step0.seq');
% rename:
% module1.mod => tipdown.mod
% module3.mod => readout.mod
% Update modules.txt accordingly
toppe.playseq(3);


return;

sys = toppe.systemspecs;

B = [1 2 3 4 42 44 48 298 300];

for ib = B
	blk = seq.getBlock(ib);
	figure;
	pulsegeq.plotblock(blk, sys.gamma, sys.raster); title(num2str(ib));
end


