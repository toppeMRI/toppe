%
%	function [g,v] = mintimegrad(Nest,g0,gf,moment,T,gmax,smax,t0,type)
%
%	
%	Function attempts to find the minimum-time gradient
%	waveform that satisfies the given design constraints.
%	Uses lpgrad -- see lpgrad for more details.
%
%	INPUT:
%		(where 	D = number of gradient dimensions,
%			Q = number of moments to rewind (0,1,2,...),
%
%	    Nest	= Starting Estimate for #points in possible solution
%				(if 2x1, Nest(1) is upper, Nest(2) is lower.	
%	    g0  	= 1xD array;  gradient value at start (G/cm).
%	    gf	   	= 1xD array;  gradient value at end (G/cm).
%	    moment 	= QxD array;  moment of gradient (/cm, s/cm, s^2/cm...)
%		
%	    T		= 1x1;  Sampling time (s).
%	    gmax	= 1x1;  Maximum gradient amplitude (G/cm).
%	    smax	= 1x1;	Maximum slew rate (G/cm/s).
%	    t0		= 1x1;  time at start of gradient (s).
%	    type	= formulation type -- 	0=lp, 
%						1=l1-norm gradient,
%						2=l1-norm gradient and slew.
%				
%	OUTPUT:
%	    g	= NxD;  Gradient waveform(s), if solution found.
% 	    v   = 1 if optimal solution found, 
%		  0 if max-iterations reached after a possible solution found
%		 -1 if no possible solution found.
%
%	Brian Hargreaves, 	Sept 2002. 

%	
%	---------------------------------------------
%	This file is maintained in CVS.
%
%	$Log: mintimegrad.m,v $
%	Revision 1.1  2018/10/25 20:39:41  jfnielse
%
%	: Committing in .
%	:
%	: Added Files:
%	: 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%	: 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%	Revision 1.1  2006/12/12 18:06:22  jfnielse
%	
%	: Added Files:
%	: 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%	: 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%	: 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%	: 	vmlpgrad.m vmsocpgrad.m
%	
%	Revision 1.7  2003/04/22 23:52:04  brian
%	Removed top/bottom check loops
%	to stream-line this code.  Added in
%	varable-cost-function for SOCP.
%	
%	Revision 1.6  2003/04/22 22:35:44  brian
%	Added comments.
%	
%	Revision 1.5  2003/04/22 18:32:55  brian
%	Added some facility for a "trisection" search
%	method, whereby the problem cost function is
%	formulated so that a solution of length N-m is
%	returned, if it exists, and the convergence
%	should then be quicker.
%	
%	Revision 1.4  2003/04/17 01:06:13  brian
%	minor edits
%	
%	Revision 1.3  2002/10/07 17:48:39  brian
%	2D spiral design w/ opt rewind.
%	
%	Revision 1.2  2002/09/17 19:21:21  bah
%	Added support for estimates for top and bottom of
%	range.  Both are checked for validity before iteration.
%	Also added arbitrary bisection point, as expected time
%	is longer to show non-existance than existance of solutions.
%	
%	Revision 1.1  2002/09/17 17:30:15  bah
%	Added to cvs
%	
%	
%	---------------------------------------------


function [g,v] = mintimegrad(Nest,g0,gf,moment,T,gmax,smax,t0,type)


n1 = clock;	% Record start time.

% ==========================================================
% 	Set Defaults (copied from lpgrad.m)
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
	type=2;
end;

if (type>2)
	method=2;
else
	method=1;	% 1=lpgrad, 2=socpgrad.
end;

maxiter=50;	% Maximum number of times that linprog will be called.
maxconslimit = 0;	% Flag - set to 1 if too many constraints.
bisectfrac = 0.6;	% Point between Nbot and Ntop for next guess.
if (type==3)
	varcost=1;	% Variable costN, number of points to zero at end.
else
  varcost=0;
  costN = 0;		% Number of points to try to zero at end.
end

Ntop=Nest(1);		% Ntop = lowest N that (ultimately) has confirmed soln.
if (length(Nest)>1)
	Nbot=Nest(2);	% Nbot = highest N that is (to be)confirmed infeasible.
else
	Nbot=.8*Ntop;	% Nbot = highest N that is confirmed infeasible.
end;

Ntopok=0;	% Set to 1 when confirmed!
Nbotok=0;	% Set to 1 when confirmed

geps = .0001;	% Test if gradient is even shorter than asked for!


% ==== Check Estimated Number of Constraints =======
if ((max(Nest) > 50) & (length(g0)>2) & (type < 2))
	v=-2;
	Nbotok=1;
	done=1;
	maxconslimit=1;
	disp('Predicted Number of Constraints is too high.');
end;

iter = 0;	% Number of iterations.

% ========================================================
% 	Newer Scheme to eliminate Top/Bottom check.
%
%	-Assume initial guesses have been checked.
%	-Start checking between Ntop and Nbot.
%	-Update Nbot and Ntop accordingly.

% =========================================================
% 	Now bisect intervals to find Nopt
%		Don't bisect evenly though, as linprog()
%		takes longer to determine that there is
%		no solution than that there is!
% =========================================================

done = 0;
if (maxconslimit==1)
	done=1;
end;

while (done==0);

	% === Set N between Nbot and Ntop. ===
	df = Ntop-Nbot;
	N = round(Nbot+df*bisectfrac);
	if (varcost==1)
		costN = floor(0.5*(Ntop-N));
	end;
	iter=iter+1;
        disp(' ');

	% === Check for solution with given N. ===
        if (method==1)
		tt = sprintf('lpgrad() - Ntop=%d  Nbot=%d   N=%d',Ntop,Nbot,N);
		disp(tt);
	        [g,v] = lpgrad(N,g0,gf,moment,T,gmax,smax,t0,type);
        else
		tt = sprintf('socpgrad() - Ntop=%d  Nbot=%d   N=%d',Ntop,Nbot,N);
		disp(tt);
	        [g,v] = socpgrad(N,g0,gf,moment,T,gmax,smax,t0,costN);
	end;


	% === If Solution Found ===
	if (v>0)
		if (Nbotok==0)	% Move bottom down if not confirmed.
			Nbot = Nbot-(Ntop-N);
			Nbot = max(Nbot,2);	% Keep at minimum of 2.
		end;
		Ntop=N;
		Ntopok=1;	% Ntop is now confirmed.
		gtop=g;		% Gradient at Ntop.

		% === Check if gradient ends earlier than N, and
		%     update Ntop and Nbot accordingly if we can. ======

		if (max(abs(gf))==0)
			Nextra = numtrailzeros(g,geps);

			if (Nextra > 1) % Extra Zeros.
				tt = sprintf('Ntop=%d, Nextra=%d, NewNtop=%d',Ntop,Nextra,Ntop+1-Nextra);
				disp(tt);
				gtop = g(1:Ntop+1-Nextra,:);	% Better gtop
				Ntop = Ntop+1-Nextra;		% Better Ntop.
				disp('>>>Extra zero grad sample(s)');
			end;

			if (Nextra < costN)  % No solution for N-costN, or
					     % we would have found it.
				Nbot = max(N-costN,Nbot);
				tt=sprintf('>>>Moving Nbot to %d',Nbot');
				disp(tt);
				Nbotok = 1;
			end;
		end;



	else   % === No Solution Found ===

		if (Ntopok==0)	% Move top up if not confirmed.
			Ntop = Ntop + (N-Nbot);
		end;
		Nbot=N;
		Nbotok=1;	% Nbot is now confirmed.
	end;
	if (Ntop<=Nbot+1)
		N=Ntop;
		g=gtop;
		done=1;
		v=1;
	end;	
	if (iter >= maxiter)
		g=gtop;
		v=0;		
		done=1;
	end;
end;


if (maxconslimit==1)
	v = -2;
	g = 0*g0;
end;

% --------- Record End time, and display times ---------

n2 = clock;	
if (1==1)
    disp(' ');
    tt=sprintf('--Start Time  %02d:%02d:%02d ',round(n1(4:6)));
    disp(tt);
    tt=sprintf('--End Time    %02d:%02d:%02d ',round(n2(4:6)));
    disp(tt);
end;
tt = sprintf('%d Points in Gradient  -- %7.4f ms',length(g),length(g)*T*1000);
disp(tt);

