% spiral readout design using Brian Hargreaves' code
%cd ~jfnielse/matlab/util/mintgrad/

mxs = 10000; % go easy on slew to avoid PNS. See also makebalanced.
mxg = 3.5;

%% DESS MWF
% Bloch-Siegert
doreverse= 0;
dovardens = 0;
fov = 22;     % cm
nleafs = 4;
npix = npix/4;

rmax = npix/(2*fov);   % max k-space radius

if dovardens
	fov0 = fov;        
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,1*(fov/20-fov)/rmax,1*(fov/40-fov)/rmax^2,rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,1000*fliplr([0.7355 -2.0987 2.1552 -0.9008 0.0869 0.0228]),rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,-1*(fov-fov/2)/rmax,-0*0.5*(fov-fov/4)/rmax^2,rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,-0.3/nleafs*fov/rmax,-0.0*fov/rmax^2,rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,-fov1/rmax,fov2/rmax^2,rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov0,-1*(fov0-fov)/rmax,-0*(fov0-fov)/rmax^2,rmax);
	%[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,-0*(fov0-fov)/rmax,-(2*fov/3)/rmax^2,rmax);
	c1 = 0.7;
	[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,-c1*fov/rmax,-(0.9-c1)*fov/rmax^2,rmax);
	cmd = sprintf('Created with Brian Hargreaves'' code: [k,g] = vds(%d,%d,4e-6,%d,fov,-%.1f*fov/rmax,-(0.9-%.1f)*fov/rmax^2,rmax);',mxs,mxg,nleafs,c1,c1);
else
	[k,g] = vds(mxs,mxg,4e-6,nleafs,fov,0,0,rmax);
	cmd = sprintf('Created with Brian Hargreaves'' code: [k,g] = vds(%d,%d,4e-6,%d,fov,0,0,rmax);',mxs,mxg,nleafs);
end
fprintf(1,'\nreadout duration is %.2f ms \n', 4e-3*size(g,2));
g = [real(g(:)) imag(g(:))];       % k = [real(k(:)) imag(k(:))]; 
g = makebalanced(g,mxs);   
g = [0 0; 0 0; g];  % add a couple of zeroes to make sure k=0 is sampled
g = makeevenlength(g);
if dovardens
	type = 'vds';
else
	type = 'unif';
end
if doreverse
	g = flipdim(g,1);  % reverse spiral 
	fname = sprintf('g-reverse-%s-nl%d-fov%d-npix%d-%s.mod', type, nleafs, fov, npix, date);
else
	fname = sprintf('g-%s-nl%d-fov%d-npix%d-%s.mod', type, nleafs, fov, npix, date);
end
fprintf(1,'gradient duration is %.2f ms \n', 4e-3*size(g,1));
rhfrsize = size(g,1);              % length of data to acquire (all of it)
addrframp = false;
addzeros = false;
desc = sprintf('spiral readout .mod file for toppe\n%s\n',cmd);
%mat2mod(0.01*ones(size(g(:,1))),0*g(:,1),g(:,1),g(:,2),0*g(:,1),90,fname,desc, 0, [0 rhfrsize],addrframp,addzeros);
toppe.writemod('gx', g(:,1), 'gy', g(:,2), 'ofname', fname, 'desc', desc);
%[desc,rho,theta,gx,gy,gz,paramsint16,paramsfloat] = readwav(fname);

return;

