%
%
%function [g,v,f,A,b,C,d,N] = 
%
%       vmsocpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0,costN)
%
%   
%   Function finds an N-point solution (if it exists) to
%   a voltage-model constrained gradient waveform design problem given 
%   second-order cone programming techniques.
%
%   This differs from vmlpgrad in that the quadratic constraints
%   on multi-dimensional gradient design can be expressed
%   as quadratic, rather than piecewise linear constraints.
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
%           Imax        = 1x1 or 1xD  -- Amplifier current limit (A)
%           Vmax        = 1x1 or 1xD  -- Amplifier voltage limit (V).
%           Rcoil       = 1x1 or 1xD  -- coil resistance in ohms.
%           Lcoil       = 1x1 or 1xD  -- coil inductance in H.
%           eta         = 1x1 or 1xD  -- coil efficiency in G/cm/A.
%       t0      = 1x1;  time at start of gradient (s).
%       type    = formulation type --   0=lp, 
%                       1=l1-norm gradient,
%                       2=l1-norm gradient and slew.
%	costN = number of points to weight at end of gradient.
%               
%   OUTPUT:
%       g   = NxD;  Gradient waveform(s), if solution found.
%       v   = 1x1;  1 if solution found, -1 if not.
%           -2 if too many constraints
%       
%       f,A,b,C,d,N = Parameters actually passed to SOCP,
%             see socp.m for more information.
%
%
%
%   Brian Hargreaves,   Oct 2002. 

%   Note 1:     This function uses SOCP.  SOCP is a development
%           from Stephen Boyd's research group at Stanford,
%           which was published around 1997.
%
%           The gradient design problem is a 
%           convex optimization problem, and can be put in
%           "standard SOCP form."  Methods, such as socp()
%           are guaranteed to find the global minimum of a
%           cost function if a feasible solution exists.
%
%
%   
%   ---------------------------------------------
%   This file is maintained in CVS.
%
%   $Log: vmsocpgrad.m,v $
%   Revision 1.1  2018/10/25 20:39:42  jfnielse
%
%   : Committing in .
%   :
%   : Added Files:
%   : 	README calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%   : 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m socpgrad.m
%   : 	tutorial.m vds.m vlim2slim.m vmlpgrad.m vmsocpgrad.m
%
%   Revision 1.1  2006/12/12 18:06:25  jfnielse
%
%   : Added Files:
%   : 	calcgradinfo.m lpgrad.m mintimegrad.m numtrailzeros.m
%   : 	plotgradinfo.m q2r21.m qdf.m slim2vlim.m socp.m
%   : 	socp_mex.mexglx socpgrad.m tutorial.m vds.m vlim2slim.m
%   : 	vmlpgrad.m vmsocpgrad.m
%
%   Revision 1.4  2003/04/22 20:16:31  brian
%   Bug fixed - elimination of redundant gradient constraints.
%
%   Revision 1.3  2003/04/22 18:26:15  brian
%   Added facility to include slack variables at the
%   end of the gradient.  These represent gradient
%   magnitude, and can be "encouraged" to go to zero
%   at the end of the gradient.
%
%   Revision 1.2  2003/04/17 17:54:48  brian
%   Changes at home to not conk out if matlab optimization toolbox is absent.
%
%   Revision 1.1  2003/04/17 01:02:04  brian
%   Added new socpgrad()
%   
%   Revision 1.4  2002/10/07 17:51:02  brian
%   try...end to take care of multi-calls
%   
%   Revision 1.3  2002/10/07 17:48:39  brian
%   2D spiral design w/ opt rewind.
%   
%   Revision 1.2  2002/10/03 01:26:16  brian
%   fixed negative g0/gf bug
%   
%   Revision 1.1  2002/10/03 00:52:04  brian
%   Added convex optimization function socp (and socp_mex.mexglx).
%   socpgrad is like lpgrad, but uses socp.
%   
%
%   Started from Rev 1.7 of lpgrad.m
%   
%   ---------------------------------------------


function [g,v,f,A,b,C,d,N] = vmsocpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0, costN)


n1 = clock; % Record start time.

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
        Imax=218;       % A.  Allow some padding.
end;
if (nargin < 7);
        Vmax = 1350;    % volts.
end;
if (nargin < 8);
       Rcoil = .35;    % ohms.
end;
if (nargin < 9);
        Lcoil = .0014;  % Henrys
end;
if (nargin < 10);
        eta = 1/56;     % G/cm/A.
end;
if (nargin < 11)
    t0=0;       % seconds.
end;
if (nargin < 12)
    costN=0;       
end;


gmax = eta*Imax;    % This is a simpler way to express the current
            % limit for equations.


% ==========================================================
%   Define D = #dimensions  and  Q = #moments
% ==========================================================

if (~isempty(g0))
    s = size(g0);
    D = s(2);
elseif (~isempty(gf))
    s = size(gf);
    D = s(2);
elseif (~isempty(moment))
    s = size(moment);
    D = s(2);
end;

%   ----- Q = # moments ------
    
if (~isempty(moment))
    s = size(moment);
    Q = s(1);
end;



% ==========================================================
%   Constants for Formulation
% ==========================================================

gamma = 4258;       % G/cm.
geps = .001;        % G/cm.  Epsilon for meeting g0 and gf constraints.
meps = [.0001;.00001;.00001;.000001;.0000001];
            % s^q/cm, Epsilon for meeting qth moment.

MaxConstraints = 600-100*D;  %

% =============================================================
%   Boundary Constraints
% =============================================================

    % ----- Initial Gradient -----

if (~isempty(g0))
    Ag0 = zeros(D,D*N);
    bg0 = zeros(D,1);
    Cg0 = zeros(D,D*N);
    dg0 = zeros(D,1);
    Ng0=ones(D,1);

    for d=1:D               % (for each dimension)
    Ag0(d,(d-1)*N+1)= 1;        % Select first point on axis...
    bg0(d,1) = -g0(1,d);        %   
    dg0(d,1)=geps;          % ...to be within geps of g0.
    end;
else
    Ag0=[]; % No constraints on initial gradient amplitude. 
    bg0=[];
    Cg0=[];
    dg0=[];
    Ng0=[];
end;


    % ----- Final Gradient -----

if (~isempty(gf))
    Agf = zeros(D,D*N);
    bgf = zeros(D,1);
    Cgf = zeros(D,D*N);
    dgf = zeros(D,1);
    Ngf=ones(D,1);

    for d=1:D               % (for each dimension)
    Agf(d,d*N)= 1;          % Select first point on axis...
    bgf(d,1) = -gf(1,d);        %   
    dgf(d,1)=geps^2;        % ...to be within geps of gf.
    end;
else
    Agf=[]; % No constraints on initial gradient amplitude. 
    bgf=[];
    Cgf=[];
    dgf=[];
    Ngf=[];
end;




    % ----- Moments -----

if (~isempty(moment))
    tm = ([1:N]-0.5)*T +t0; % Array of center times for gradient points.

    tt = [ones(size(tm)); tm; tm.^2+T^2/12; tm.^3+T^2/2*tm]*gamma*T;  

                % tt are multiplied by g to discretely
                % approximate the gradient moments.

    Am = zeros(2*D*Q,D*N);
    bm = zeros(2*D*Q,1);
    Cm = zeros(2*D*Q,D*N);
    dm = zeros(2*D*Q,1);
    Nm = ones(2*D*Q,1);

    if (Q>4)
    disp('4th moments and higher are not supported');
    end;
    for d=1:D
        for q=1:Q   
            Cm(2*d-1+(q-1)*2*D ,(d-1)*N+1:d*N) = -tt(q,:); %negate wrt lpgrad.
            Cm(2*d  +(q-1)*2*D ,(d-1)*N+1:d*N) =  tt(q,:);
            dm(2*d-1+(q-1)*2*D ,1) =  moment(q,d)+meps(q);
            dm(2*d  +(q-1)*2*D ,1) = -moment(q,d)+meps(q);
    end;
    end;
else
    Am=[];  % No constraints on moments.
    bm=[];
    Cm=[];  
    dm=[];
    Nm = [];
end;



% =============================================================
%   Gradient Amplitude Constraints
% =============================================================
%
%   Note 3: We must generally constrain each gradient
%       point so that the gradient vector amplitude is
%       less than gmax.  For SOCP, we use quadratic constraints
%       for any number of dimensions.
%
%       Furthermore, if the endpoints are constrained, then
%       the slew rate constraint means that for some number of
%       points at each end, some gradient amplitude constraints
%       are unnecessary.
%



%   ------  Some constraints are redundant based on max
%       slew constraints and boundary constraints -------

igstart = 0;    % # points at start where we can ignore gmax constraint.
igend = 0;  % # points at end...
    

if (1==1)   % Eliminate Redundant Gradient constraints.
        % These are constraints that could not be violated
        % since the boundary point is too far the constraint to be
        % reached in the given # of points.


    smax = (eta*Vmax+Rcoil*gmax)/Lcoil; % Max POSSIBLE slew rate.

    if (~isempty(g0))
        g0dist = gmax-norm(g0);
        igstart = floor(abs(g0dist/T/smax));    
    end;
    if (~isempty(gf))
        gfdist = gmax-norm(gf);
        igend = floor(abs(gfdist/T/smax));  
    end;
end;

% -- disable elimination of constraints. ---
if (1==0)
    igstart=0;
    igend=0;
end;

igend = N+1-igend;  % Convert to point number.

    
%   ------ Now expand the constraints to apply to each point
%   
ngcons=0;
Agcons=zeros(N*D,N*D);

for k=1:N
    if ((k >= igstart) & (k <= igend))  % be conservative!
        ngcons=ngcons+1;
        for d=1:D
            Agcons((ngcons-1)*D+d,(d-1)*N+k)=1; 
        end;
    end;
end;
Agcons = Agcons(1:D*ngcons,:);  % Discard zero rows at end.
bgcons = zeros(D*ngcons,1);
Cgcons = zeros(ngcons,N*D); 
dgcons = gmax*ones(ngcons,1);
Ngcons = D*ones(ngcons,1);      



% =============================================================
%   Slew Constraints
% =============================================================
%
%   Note 4: These are similar to the gradient constraints (Note 3).
%       However, there are only N-1 "points", as the DIFFERENCE
%       between gradient points is what is constrained.
%
%       Similar to gradient constraints, these are |g[k]-g[k-1]|^2
%       for any # dimensions.
%
 



Ascons = zeros((N-1)*D, N*D);
bscons = zeros((N-1)*D, 1);
Cscons = zeros((N-1), N*D);
%dscons = smax*T* ones((N-1), 1);   OLD CSRL Constraint.
dscons = Vmax*T/Lcoil* ones((N-1), 1);
Nscons =  D*ones((N-1), 1);

alpha = Rcoil*T/2/Lcoil-1;  % Voltage-Limit constraint for use below.
beta = Rcoil*T/2/Lcoil+1;   % Voltage-Limit constraint for use below.

for k=1:(N-1)
    for d=1:D
        Ascons((k-1)*D+d,(d-1)*N+k   )=    alpha;
        Ascons((k-1)*D+d,(d-1)*N+k+1 )=    beta;
    end;
end;

% ============================================================
%   Cost Slack Constraints
% ============================================================
%
%   Note 4.5:  Here we add some slack variables so that
%	we can use a cost function that tries to shorten
%	the number of non-zero gradient points.  Each slack
%	variable is H >= sqrt(Gx^2+Gy^2).  

NNslack = costN;
if (NNslack > 0)
	NNslack = min(NNslack,N);
	Aslack = zeros(D*NNslack, N*D+NNslack);	% Will update.
	bslack = zeros(D*NNslack, 1);
	Cslack = zeros(NNslack, N*D+NNslack);   % Will update.
	dslack = zeros(NNslack, 1);
	Nslack = D*ones(NNslack,1);

	for k=1:NNslack
	    for d=1:D
		Aslack((k-1)*D+d,N*d-NNslack+k) = 1;
	    end;
	    Cslack(k,N*D+k) = 1;
	end;
end;

% ============================================================
%       Combine All Constraints
% ============================================================
%
%   Note 5: Here it is worth looking back at Notes 1 and 2.
%   The constraints are all put together to constrain the
%   X variable.



    if (NNslack == 0)
    	AA = [Ag0; Agf; Am; Agcons; Ascons];
    	bb = [bg0; bgf; bm; bgcons; bscons];
    	CC = [Cg0; Cgf; Cm; Cgcons; Cscons];
    	dd = [dg0; dgf; dm; dgcons; dscons];
    	NN = [Ng0; Ngf; Nm; Ngcons; Nscons];
    	s = size(AA); 
    	ff = zeros(s(2),1);
    	ff(1:N)=1;
    else
    	AA = [Ag0; Agf; Am; Agcons; Ascons];
	sAA = size(AA);
	%size(Aslack)
	AA = [[AA zeros(sAA(1),NNslack)]; Aslack];
    	bb = [bg0; bgf; bm; bgcons; bscons; bslack];
    	CC = [Cg0; Cgf; Cm; Cgcons; Cscons];
	sCC = size(CC);
	CC = [[CC zeros(sCC(1),NNslack)]; Cslack];
    	dd = [dg0; dgf; dm; dgcons; dscons; dslack];
    	NN = [Ng0; Ngf; Nm; Ngcons; Nscons; Nslack];
	ff = [zeros(N*D,1); ones(NNslack,1)];

    end;
    

s=size(AA);
tt = sprintf('Calling socp() with %d variables and %d constraints',s(2),length(NN));
disp(tt);

if (1==0)
    AA
    bb
    CC
    dd
    NN
    N
    rank([AA;CC])
    rank(AA)
end;

% ----------- Call socp() ------------------------
%
%  Note:    socp() generates an error if there is no solution.
%       Thus we use the try...catch...end formulation!
%
if (length(NN) <= MaxConstraints)
  if (exist('socp_mex')~=0)
    
    try
      v=1;
      [g,info] = socp(ff,AA,bb,CC,dd,NN);
    catch
      v=-1;
      disp('No Solution');
    end;
  else
    sz = size(AA);
    g = zeros(sz(2),1);
    v = 1;
    disp('*** socp_mex does not exist.  Returning zero solution. ***');
  end;
else
    g = 0*g0;
    v = -2;
end;


if (v>0)
    g = reshape(g(1:N*D),N,D);
else
    g = 0*g0;
end;

% --------- Record End time, and display times ---------

n2 = clock; 
if (0==1)
    tt=sprintf('Start Time  %02d:%02d:%02d ',round(n1(4:6)));
    disp(tt);
    tt=sprintf('End Time    %02d:%02d:%02d ',round(n2(4:6)));
    disp(tt);
end;

% --------- Set output arguments ----------


f=ff;
A=AA;
b=bb;
C=CC;
d=dd;
N=NN;

