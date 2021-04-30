function plotblock(blk, gamma, raster)
% 3T: gamma = 4257.6 Hz/G
% raster: 4e-6; 

rf = [];
gx = [];
gy = [];
gz = [];

if ~isempty(blk.gx)
	gx = pulsegeq.sub_trap2shape(blk.gx, gamma, raster);
end
if ~isempty(blk.gy)
	gy = pulsegeq.sub_trap2shape(blk.gy, gamma, raster);
end
if ~isempty(blk.gz)
	gz = pulsegeq.sub_trap2shape(blk.gz, gamma, raster);
end
	
plot([gx gy gz], 'o');
ylabel('G/cm');
legend('gx', 'gy', 'gz');
