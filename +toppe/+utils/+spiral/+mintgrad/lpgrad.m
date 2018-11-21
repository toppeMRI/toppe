%
%   function g = lpgrad(N,g0,gf,moment,T,gmax,smax,t0,type)
%
%   
%   Function does the same as vmlpgrad, but uses a 
%   constant-slew-rate-limit model.
%
%   lpgrad uses gmax and smax.
%   ------
%   vmlpgrad uses Imax, Vmax, Rcoil, Lcoil, eta 
%
%   
%   Function finds an N-point solution (if it exists) to
%   a constrained gradient waveform design problem given 
%   linear programming techniques.
%
%   INPUT:
%       (where  D = number of gradient dimensions,
%           Q = number of moments to rewind (0,1,2,...),
%
%       N       = Number of points in possible solution.    
%       g0      = 1xD array;  gradient value at start (G/cm).
%       gf      = 1xD array;  gradient value at end (G/cm).
%       moment  = QxD array;  moment of gradient (/cm, s/cm, s^2/cm...)
%       
%       T       = 1x1;  Sampling time (s).
%       gmax    = 1x1;  Maximum gradient amplitude (G/cm).
%       smax    = 1x1;  Maximum slew rate (G/cm/s).
%       t0      = 1x1;  time at start of gradient (s).
%       type    = formulation type --   0=lp, 
%                       1=l1-norm gradient,
%                       2=l1-norm gradient and slew.
%               
%   OUTPUT:
%       g   = NxD;  Gradient waveform(s), if solution found.
%       v   = 1x1;     Value returned from linprog().  -1,0 or 1.
%       f   =   Cost function passed to linprog().
%       A   =   Constraint matrix passed to linprog().
%       b   =   Constraint vector passed to linprog().
% 
%
%   Brian Hargreaves,   Sept 2002. 

%   Note 1:     This function uses linear programming functions
%           such as linprog() that solve linear optimization
%           problems.  The gradient design problem is a 
%           convex optimization problem, and can be put in
%           "standard LP form."  Methods, such as linprog()
%           are guaranteed to find the global minimum of a
%           cost function if a feasible solution exists.
%
%
%   
%   ---------------------------------------------
%   This file is maintained in CVS.
%
%   $Log: lpgrad.m,v $
%   Revision 1.1  2018/10/25 20:39:41  jfnielse
%
%   : Committing in .
%   :
%   : Added Files:
%   : 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%   : 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%   : 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%   Revision 1.1  2006/12/12 18:06:21  jfnielse
%
%   : Added Files:
%   : 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%   : 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%   : 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%   : 	vmlpgrad.m vmsocpgrad.m
%
%   Revision 1.11  2003/04/17 17:54:48  brian
%   Changes at home to not conk out if matlab optimization toolbox is absent.
%
%   Revision 1.10  2003/04/17 01:02:57  brian
%   Modified so that these use vmsocpgrad.m and vmlpgrad.m,
%   the voltage models.  The goal is to reduce the amount
%   of code that is maintained!!
%   
%   Revision 1.9  2002/10/07 17:49:18  brian
%   Added error checking on g0/gf
%   
%   Revision 1.8  2002/10/07 17:48:39  brian
%   2D spiral design w/ opt rewind.
%   
%   Revision 1.7  2002/10/02 22:19:04  brian
%   --Made tolerances smaller.
%   --Fixed A1 bug
%   --Added code to try to force shorter waveform.
%   
%   Revision 1.6  2002/09/19 18:24:34  bah
%   Added support for 3D designs with type 2 formulation.
%   
%   Revision 1.5  2002/09/17 17:29:37  bah
%   Added numerous comments
%   
%   Revision 1.4  2002/09/17 02:25:24  bah
%   l1-norm slew constraint added
%   
%   Revision 1.3  2002/09/16 21:45:56  bah
%   Fixed moment-constraint bug so axes are indepently constrained
%   
%   Revision 1.2  2002/09/16 19:26:27  bah
%   Added l1 formulation for gradient.
%   Still testing... not so successfully yet!
%   
%   Revision 1.1  2002/09/16 16:30:51  bah
%   added to cvs.
%   
%   ---------------------------------------------


function [g,v,f,A,b] = lpgrad(N,g0,gf,moment,T,gmax,smax,t0,type)


% ==========================================================
%   Set Defaults
% ==========================================================

if (nargin < 4)     
    moment=[];
end;
if (nargin < 5);
    T=0.000004; % seconds.
end;
if (nargin < 6);
    gmax=3.9;   % G/cm.  Allow some padding.
end;
if (nargin < 7);
    smax=14500; % G/cm/s.  Allow some padding for quantization.
end;
if (nargin < 8)
    t0=0;       % seconds.
end;
if (nargin < 9)
    type=2;
end;

[Imax,Vmax,Rcoil,Lcoil,eta] = slim2vlim(gmax,smax);
[g,v,f,A,b] = vmlpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0,type);
