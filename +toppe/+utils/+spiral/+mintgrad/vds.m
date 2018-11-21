%
%	function [k,g,s,time,r,theta] = vds(smax,gmax,T,N,F0,F1,F2,rmax)
%
%
%	VARIABLE DENSITY SPIRAL GENERATION:
%	----------------------------------
%
%	Function generates variable density spiral which traces
%	out the trajectory
%				 
%			k(t) = r(t) exp(i*q(t)), 		[1]
%
%	Where q is the same as theta...
%		r and q are chosen to satisfy:
%
%		1) Maximum gradient amplitudes and slew rates.
%		2) Maximum gradient due to FOV, where FOV can
%		   vary with k-space radius r, as
%
%			FOV(r) = F0 + F1*r + F2*r*r 		[2]
%
%
%	INPUTS:
%	-------
%	smax = maximum slew rate G/cm/s
%	gmax = maximum gradient G/cm (limited by Gmax or FOV)
%	T = sampling period (s) for gradient AND acquisition.
%	N = number of interleaves.
%	F0,F1,F2 = FOV coefficients with respect to r - see above.
%	rmax= value of k-space radius at which to stop (cm^-1).
%		rmax = 1/(2*resolution);
%
%
%	OUTPUTS:
%	--------
%	k = k-space trajectory (kx+iky) in cm-1.
%	g = gradient waveform (Gx+iGy) in G/cm.
%	s = derivative of g (Sx+iSy) in G/cm/s.
%	time = time points corresponding to above (s).
%	r = k-space radius vs time (used to design spiral)
%	theta = atan2(ky,kx) = k-space angle vs time.
%
%
%	METHODS:
%	--------
%	Let r1 and r2 be the first derivatives of r in [1].	
%	Let q1 and q2 be the first derivatives of theta in [1].	
%	Also, r0 = r, and q0 = theta - sometimes both are used.
%	F = F(r) defined by F0,F1,F2.
%
%	Differentiating [1], we can get G = a(r0,r1,q0,q1,F)	
%	and differentiating again, we get S = b(r0,r1,r2,q0,q1,q2,F)
%
%	(functions a() and b() are reasonably easy to obtain.)
%
%	FOV limits put a constraint between r and q:
%
%		dr/dq = N/(2*pi*F)				[3]	
%
%	We can use [3] and the chain rule to give 
%
%		q1 = 2*pi*F/N * r1				[4]
%
%	and
%
%		q2 = 2*pi/N*dF/dr*r1^2 + 2*pi*F/N*r2		[5]
%
%
%
%	Now using [4] and [5], we can substitute for q1 and q2
%	in functions a() and b(), giving
%
%		G = c(r0,r1,F)
%	and 	S = d(r0,r1,r2,F,dF/dr)
%
%
%	Using the fact that the spiral should be either limited
%	by amplitude (Gradient or FOV limit) or slew rate, we can
%	solve 
%		|c(r0,r1,F)| = |Gmax|  				[6]
%
%	analytically for r1, or
%	
%	  	|d(r0,r1,r2,F,dF/dr)| = |Smax|	 		[7]
%
%	analytically for r2.
%
%	[7] is a quadratic equation in r2.  The smaller of the 
%	roots is taken, and the real part of the root is used to
%	avoid possible numeric errors - the roots should be real
%	always.
%
%	The choice of whether or not to use [6] or [7], and the
%	solving for r2 or r1 is done by q2r21.m.
%
%	Once the second derivative of theta(q) or r is obtained,
%	it can be integrated to give q1 and r1, and then integrated
%	again to give q and r.  The gradient waveforms follow from
%	q and r. 	
%
%	Brian Hargreaves -- Sept 2000.
%
%	See Brian's journal, Vol 6, P.24.
%
%
%	See also:  vds2.m
%

% =============== CVS Log Messages ==========================
%	$Log: vds.m,v $
%	Revision 1.3  2018/10/26 01:38:37  jfnielse
%
%	: Committing in .
%	:
%	: Modified Files:
%	: 	makesosreadout.m +mintgrad/vds.m
%
%	Revision 1.2  2018/10/26 01:12:18  jfnielse
%
%	: Committing in .
%	:
%	: Modified Files:
%	: 	q2r21.m vds.m
%
%	Revision 1.1  2018/10/25 20:39:41  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%	: 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%	Revision 1.1  2006/12/12 18:06:24  jfnielse
%	
%	: Added Files:
%	: 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%	: 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%	: 	vmlpgrad.m vmsocpgrad.m
%	
%	Revision 1.3  2002/11/18 05:36:02  brian
%	Rounds lengths to a multiple of 4 to avoid
%	frame size issues later on.
%	
%	Revision 1.2  2002/11/18 05:32:19  brian
%	minor edits
%	
%	Revision 1.1  2002/03/28 01:03:20  bah
%	Added to CVS
%	
%
% ===========================================================

function [k,g,s,time,r,theta] = vds(smax,gmax,T,N,F0,F1,F2,rmax)

import toppe.utils.spiral.mintgrad.*

%disp('vds.m');
gamma = 4257.6;

oversamp = 8;		% Keep this even.
To = T/oversamp;	% To is the period with oversampling.



q0 = 0;	
q1 = 0;
theta = zeros(1,10000);
r = zeros(1,10000);
r0 = 0;
r1 = 0;

time = zeros(1,10000);
t = 0;
count = 1;

theta = zeros(1,1000000);
r = zeros(1,1000000);
time = zeros(1,1000000);

while r0 < rmax
	[q2,r2] = q2r21(smax,gmax,r0,r1,To,T,N,[F0,F1,F2]);

	% Integrate for r, r', theta and theta' 	
	q1 = q1 + q2*To;
	q0 = q0 + q1*To;
 	t = t + To;

	r1 = r1 + r2*To;
	r0 = r0 + r1*To;

	% Store.
	count = count+1; 
	theta(count) = q0;
	r(count) = r0;
	time(count) = t;

	if (rem(count,100)==0)
		tt = sprintf('%d points, |k|=%f',count,r0);
		%disp(tt);
	end;
end;

r = r(oversamp/2:oversamp:count);
theta = theta(oversamp/2:oversamp:count);
time = time(oversamp/2:oversamp:count);

%	Keep the length a multiple of 4, to save pain...!
%
ltheta = 4*floor(length(theta)/4);
r=r(1:ltheta);
theta=theta(1:ltheta);
time=time(1:ltheta);

%
% 	Plot.
%
%x = alpha*theta .* cos(theta);
%y = alpha*theta .* sin(theta);

%plot(x,y);
%title('k-space trajectory.');


k = r.*exp(i*theta);

g = 1/gamma*([k 0]-[0 k])/T;
g = g(1:length(k));

s = ([g 0]-[0 g])/T;
s = s(1:length(k));


% ========= Plot gradients and slew rates. ==========


subplot(2,2,1);
plot(real(k),imag(k));
title('k_y vs k_x');
axis('square');

subplot(2,2,2);
plot(time,real(k),'r--',time,imag(k),'b--',time,abs(k),'k-');
title('k vs t');
ylabel('k (cm^{-1})');

subplot(2,2,3);
plot(time,real(g),'r--',time,imag(g),'b--',time,abs(g),'k-');
title('g vs t');
ylabel('G (G/cm)');

subplot(2,2,4);
plot(time,real(s),'r--',time,imag(s),'b--',time,abs(s),'k-');
title('s vs t');
ylabel('Slew Rate (G/cm/s)');

