%
%
%	function [g,v,f,A,b,C,d,N] 
%	
%			= socpgrad(N,g0,gf,moment,T,gmax,smax,t0,costN)
%
%
%	Function finds an N-point solution (if it exists) to
%	a constrained gradient waveform design problem given 
%	second-order cone programming techniques.
%
%	This differs from lpgrad in that the quadratic constraints
%	on multi-dimensional gradient design can be expressed
%	as quadratic, rather than piecewise linear constraints.
%
%	INPUT:
%		(where 	D = number of gradient dimensions,
%			Q = number of moments to rewind (0,1,2,...),
%
%	    N		= Number of points in possible solution.	
%	    g0  	= 1xD array;  gradient value at start (G/cm).
%	    gf	   	= 1xD array;  gradient value at end (G/cm).
%	    moment 	= QxD array;  moment of gradient (/cm, s/cm, s^2/cm...)
%		
%	    T		= 1x1;  Sampling time (s).
%	    gmax	= 1x1;  Maximum gradient amplitude (G/cm).
%	    smax	= 1x1;	Maximum slew rate (G/cm/s).
%	    t0		= 1x1;  time at start of gradient (s).
%	    costN	= number of slack variables to add to try to
%			   shorten gradient.
%	OUTPUT:
%	    g	= NxD;  Gradient waveform(s), if solution found.
%	    v   = 1x1;  1 if solution found, -1 if not.
%			-2 if too many constraints
% 
%	    f,A,b,C,d,N == parameters actually passed to SOCP.
%
%	Brian Hargreaves, 	Oct 2002. 

%	Note 1:		This function uses SOCP.  SOCP is a development
%			from Stephen Boyd's research group at Stanford,
%			which was published around 1997.
%
%			The gradient design problem is a 
%			convex optimization problem, and can be put in
%			"standard SOCP form."  Methods, such as socp()
%			are guaranteed to find the global minimum of a
%			cost function if a feasible solution exists.
%
%
%	
%	---------------------------------------------
%	This file is maintained in CVS.
%
%	$Log: socpgrad.m,v $
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
%	Revision 1.6  2003/04/22 18:24:30  brian
%	minor edits
%	
%	Revision 1.5  2003/04/17 01:02:57  brian
%	Modified so that these use vmsocpgrad.m and vmlpgrad.m,
%	the voltage models.  The goal is to reduce the amount
%	of code that is maintained!!
%	
%	Revision 1.4  2002/10/07 17:51:02  brian
%	try...end to take care of multi-calls
%	
%	Revision 1.3  2002/10/07 17:48:39  brian
%	2D spiral design w/ opt rewind.
%	
%	Revision 1.2  2002/10/03 01:26:16  brian
%	fixed negative g0/gf bug
%	
%	Revision 1.1  2002/10/03 00:52:04  brian
%	Added convex optimization function socp (and socp_mex.mexglx).
%	socpgrad is like lpgrad, but uses socp.
%	
%
%	Started from Rev 1.7 of lpgrad.m
%	
%	---------------------------------------------


function [g,v,f,A,b,C,d,N] = socpgrad(N,g0,gf,moment,T,gmax,smax,t0,costN)


% ==========================================================
% 	Set Defaults
% ==========================================================

if (nargin < 4)		
	moment=[];
end;
if (nargin < 5);
	T=0.000004;	% seconds.
end;
if (nargin < 6);
	gmax=3.9;	% G/cm.  Allow some padding.
end;
if (nargin < 7);
	smax=14500;	% G/cm/s.  Allow some padding for quantization.
end;
if (nargin < 8)
	t0=0;		% seconds.
end;
if (nargin < 9)
	costN=0;		
end;


% Use the voltage-model version here, setting some
% parameters so that it looks like a constant-slew-rate-limit model.
%

[Imax,Vmax,Rcoil,Lcoil,eta] = slim2vlim(gmax,smax);
[g,v,f,A,b,C,d,N] = vmsocpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0,costN);


