files = dir('*.mod');
for ii = 1:length(files)
	toppe.plotmod(files(ii).name);
end
