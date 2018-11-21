%
%	function [k,g,s,m1,m2,t,v]=calcgradinfo(g,T,k0,R,L,eta)
%
%	Function calculates gradient information
%
%	INPUT:
%		g	gradient (G/cm) = gx + i*gy
%		T	sample period (s)
%		k0	initial condition for k-space.
%		R	coil resistance (ohms)
%		L	coil inductance (H)
%		eta	coil efficiency (G/cm/A)
%
%	OUTPUT:
%		k	k-space trajectory (cm^(-1))
%		g	gradient (G/cm)
%		s	slew rate trajectory (G/cm/s)
%		m1	first moment trajectory (s/cm)
%		m2	second moment trajectory (s^2/cm)
%		t	vector of time points (s)
%		v	voltage across coil.
%

%	B.Hargreaves, Aug 2002.
%

% =============== CVS Log Messages ==========================
%	$Log: calcgradinfo.m,v $
%	Revision 1.1  2018/10/25 20:39:41  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%	: 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%	Revision 1.1  2006/12/12 18:06:21  jfnielse
%	
%	: Added Files:
%	: 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%	: 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%	: 	vmlpgrad.m vmsocpgrad.m
%	
%	Revision 1.4  2002/10/14 15:59:36  brian
%	Added calculation/plot of voltage.
%	
%	Revision 1.3  2002/09/26 21:16:39  bah
%	Fixed m2 calculation
%	
%	Revision 1.2  2002/09/16 22:35:51  bah
%	now supports Nx2 g
%	
%	Revision 1.1  2002/09/05 18:35:43  bah
%	Separted into calculation and plot parts.
%		
%
% ===========================================================


function [k,g,s,m1,m2,t,v]=calcgradinfo(g,T,k0,R,L,eta)


if (nargin < 2)
	T = .000004;
end;
if (nargin < 3)
	k0 = 0;
end;
if (nargin < 4)
	R = .5;	 	% Rough estimate.
end;
if (nargin < 5)
	L = .001;	% (Estimated)
end;
if (nargin < 6)
	eta = 1/50;	% (Estimate)
end;

gamma = 4258;

s = size(g);
lg = max(size(g));
k = k0 + cumsum(g)*gamma*T;
t = ([1:length(g)]-.5)*T;
t = t';
s = size(g);
tt = t*ones(1,s(2));
s = [g; g(lg,:)]-[0*g(1,:); g];
sz=size(s);
s = s(2:sz(1),:)/T;
m1 = cumsum(g.*tt)*gamma*T;
m2 = cumsum(g.*(tt.*tt+T^2/12))*gamma*T;
v = (1/eta)*(L*s+R*g);






