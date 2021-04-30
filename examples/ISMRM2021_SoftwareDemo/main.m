%setup;

seq01_Fid;
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly

seq02_Fid_multipleFAs;
% rename:
% module1.mod => tipdown.mod
% module2.mod => readout.mod
% Update modules.txt accordingly


%seq05_RadialSE;
pulsegeq.seq2ge('se_radial.seq');
% rename:
% module1.mod => tipdown.mod
% module6.mod => readout.mod
% Update modules.txt accordingly
toppe.playseq(4,'drawpause',false);

return;

sys = toppe.systemspecs;

B = [1 2 3 4 42 44 48 298 300];

for ib = B
	blk = seq.getBlock(ib);
	figure;
	pulsegeq.plotblock(blk, sys.gamma, sys.raster); title(num2str(ib));
end


