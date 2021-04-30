%setup;

%seq05_RadialSE;
pulsegeq.seq2ge('se_radial.seq');
toppe.playseq(4,'drawpause',false);

return;

sys = toppe.systemspecs;

B = [1 2 3 4 42 44 48 298 300];

for ib = B
	blk = seq.getBlock(ib);
	figure;
	pulsegeq.plotblock(blk, sys.gamma, sys.raster); title(num2str(ib));
end


