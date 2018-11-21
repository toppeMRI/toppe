%	function [r1,r2] = qdf(a,b,c)
%
%	Outputs quadratic roots of ax^2+bx+c = 0.
%


% =============== CVS Log Messages ==========================
%	This file is maintained in CVS version control.
%
%	$Log: qdf.m,v $
%	Revision 1.1  2018/10/25 20:39:41  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%	: 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%	Revision 1.1  2006/12/12 18:06:23  jfnielse
%	
%	: Added Files:
%	: 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%	: 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%	: 	vmlpgrad.m vmsocpgrad.m
%	
%	Revision 1.1  2002/03/28 01:27:46  bah
%	Added to CVS
%	
%
% ===========================================================


function [roots] = qdf(a,b,c)

d = b^2 - 4*a*c;

roots(1) = (-b + sqrt(d))/(2*a);
roots(2) = (-b - sqrt(d))/(2*a);




