pfile = 'P,gre.7';

d = toppe.utils.loadpfile(pfile);  % [nx ncoils nz 1 ny]

d = permute(d, [1 5 3 2 4]);   % [nx ny nz ncoils]

% size(ims) = [nx ny nz ncoils]
% size(rss) = [nx ny nz]  (root sum of squares coil-combined image)
[ims, rss] = toppe.utils.ift3(d);

im(rss);  % requires MIRT toolbox (Jeff Fessler)
