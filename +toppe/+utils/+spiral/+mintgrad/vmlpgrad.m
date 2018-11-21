%
%   function g = vmlpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0,type)
%
%   
%   Function finds an N-point solution (if it exists) to
%   a constrained gradient waveform design problem given 
%   linear programming techniques.
%
%   Same as lpgrad, but uses a voltage model instead of the
%   constant Smax model.
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
%       Imax    = 1x1 or 1xD  -- Amplifier current limit (A)
%       Vmax    = 1x1 or 1xD  -- Amplifier voltage limit (V).
%       Rcoil   = 1x1 or 1xD  -- coil resistance in ohms.
%       Lcoil   = 1x1 or 1xD  -- coil inductance in H.
%       eta     = 1x1 or 1xD  -- coil efficiency in G/cm/A.
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
%   $Log: vmlpgrad.m,v $
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
%   Revision 1.6  2003/04/17 17:54:48  brian
%   Changes at home to not conk out if matlab optimization toolbox is absent.
%
%   Revision 1.5  2003/04/17 01:05:52  brian
%   minor edits
%   
%   Revision 1.4  2002/10/12 00:02:49  brian
%   Bug fixes.
%   
%   Revision 1.3  2002/10/11 22:32:49  brian
%   ready to test...
%   
%   Revision 1.2  2002/10/11 22:06:14  brian
%   --Added in defaults for new parameters.
%   --Modified elimination of constraints for gmax checks.
%   
%   Revision 1.1  2002/10/11 18:01:45  brian
%   Started to change header
%   
%
%   ------------------------------------
%   Copied from this version of lpgrad.m
%   ------------------------------------
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


function [g,v,f,A,b] = vmlpgrad(N,g0,gf,moment,T,Imax,Vmax,Rcoil,Lcoil,eta,t0,type)


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
    Imax=218;   % A.  Allow some padding.
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
    eta = 1/56; % G/cm/A.
end;
if (nargin < 11)
    t0=0;       % seconds.
end;
if (nargin < 12)
    type=2;
end;

gmax = Imax*eta;    % Still seem to use gmax a lot here...

% ==========================================================
%   Linear-Programming Formulation Type
% ==========================================================
%
%   Note 2: The basic formulation (type 0) will have all 
%       constraints expressed as linear inequalities
%       that are linear functions of the gradient waveform.
%       In this case, the cost function will often be irrelevant,
%       as it is difficult to have a useful cost function that is
%       linear, so the problem is simply one of existance.  In this
%       case the solution vector is X = [gx;gy;gz] where gx, gy and
%       gz are Nx1 gradient waveform timepoints for each axis (if 
%       2D or 3D.
%
%       An alternative (type 1) is to augment the solution vector as
%       X = [gx; gy; gz; hx; hy; hz] where h=|g| in the solution.
%       Each h is constrained so that h >= |g| by constraints
%       g-h<0 and -g-h<0, and the sum of h's is minimized in the
%       cost function.  This can be advantageous as minimizing
%       sum|g| is generally a good thing.  Also, later points
%       can be weighted more heavily in the cost function to
%       "encourage" the non-zero gradient length to be minimum.
%
%       A second alternative (type 2) is to augment the solution vector 
%       further, as X = [gx; gy; gz; hx; hy; hz; sx; sy; sz]
%       where each  s=|g[k]-g[k-1]| in the solution.  Thus there
%       are N-1 s variables for each dimension.  Now the sum
%       of h's + s's is minimized. 
%
%       Type 1 and 2 formulations have more variables in the 
%       solution, but for multi-dimensional problems can reduce
%       the number of constraints significantly, so are desirable. 
%

if (type==0)
    l1grad=0;
    l1slew=0;
elseif (type==1)
    l1grad = 1;
    l1slew = 0;
elseif (type==2)
    l1grad = 1;
    l1slew = 1;
end;


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
npwl = 4;       % # piecewise linear segments per quadrant.


% =============================================================
%   Check Boundary Constraints are within Gmax
% =============================================================

%   It is possible (and not uncommon) for the boundary
%   constraints on g to violate the gmax constraint.  This
%   should cause an error condition, as otherwise there
%   will never be a feasible solution.

err0=0;
errf=0;
if (D==1)
    if ((~isempty(g0)) & (abs(g0)>=gmax))
        err0 = 1;
    end;
    if (abs(gf)>=gmax)
        errf = 1;
    end;

else 
    if (D==2)
    [Ag,bg] = pwlapprox(gmax, npwl,1);
    elseif (D==3)
    [Ag,bg] = pwl3d9approx(gmax);
    end;    
    if (~isempty(g0))
        err0= 1;
        if (Ag*g0.'<bg)
            err0 = 0;
        end;
    end;
    if (~isempty(gf))
        errf= 1;
        if (Ag*gf.'<bg)
            errf = 0;
        end;
    end;
end;
if ((err0==1) & (errf==1))
    error('Both initial and final gradient constraints violate maximum.');
end;
if (err0==1)
    error('Initial gradient constraint violates maximum gradient constraint.');
end;
if (errf==1)
    error('Final gradient constraint violates maximum gradient constraint.');
end;





% =============================================================
%   Boundary Constraints
% =============================================================

    % ----- Initial Gradient -----

if (~isempty(g0))
    Ag0 = zeros(2*D,D*N);
    bg0 = zeros(2*D,1);
    for d=1:D               % (for each dimension)
    Ag0(2*d-1,(d-1)*N+1)= 1;    % Select first point on axis...
    Ag0(2*d  ,(d-1)*N+1)=-1;    %   "             "
    bg0(2*d-1,1) =  g0(1,d)+geps;   % ...to be within geps of g0.
    bg0(2*d  ,1) = -g0(1,d)+geps;   %   "         "
    end;
else
    Ag0=[]; % No constraints on initial gradient amplitude. 
    bg0=[];
end;


    % ----- Final Gradient -----

if (~isempty(gf))
    Agf = zeros(2*D,D*N);
    bgf = gf.';
    for d=1:D               % See comments above for initial
    Agf(2*d-1,d*N)= 1;      % gradient value constraint.
    Agf(2*d  ,d*N)=-1;      %  (this is the same...)    
    bgf(2*d-1,1) =  gf(1,d)+geps;   
    bgf(2*d  ,1) = -gf(1,d)+geps;
    end;
else
    Agf=[]; % No constraints on final gradient. 
    bgf=[];
end;


    % ----- Moments -----

if (~isempty(moment))
    tm = ([1:N]-0.5)*T +t0; % Array of center times for gradient points.

    tt = [ones(size(tm)); tm; tm.^2+T^2/12; tm.^3+T^2/2*tm]*gamma*T;  

                % tt are multiplied by g to discretely
                % approximate the gradient moments.
    Am = zeros(2*D*Q,D*N);
    bm = zeros(2*D*Q,1);
    if (Q>4)
    disp('4th moments and higher are not supported');
    end;
    for d=1:D
        for q=1:Q   
            Am(2*d-1+(q-1)*2*D ,(d-1)*N+1:d*N) =  tt(q,:);
            Am(2*d  +(q-1)*2*D ,(d-1)*N+1:d*N) = -tt(q,:);
            bm(2*d-1+(q-1)*2*D ,1) =  moment(q,d)+meps(q);
            bm(2*d  +(q-1)*2*D ,1) = -moment(q,d)+meps(q);
    end;
    end;
else
    Am=[];  % No constraints on moments.
    bm=[];
end;



% =============================================================
%   Gradient Amplitude Constraints
% =============================================================
%
%   Note 3: We must generally constrain each gradient
%       point so that the gradient vector amplitude is
%       less than gmax.  In 1D, this is just the magnitude.
%       In 2D or 3D, we approximate this (quadratic) constraint
%       by many linear constraints.
%
%       If the formulation type is 1 or 2, then |g| is
%       available, and the number of constraints can be reduced
%       to those where the constraint coefficients on g's are
%       positive.
%
%       Furthermore, if the endpoints are constrained, then
%       the slew rate constraint means that for some number of
%       points at each end, some gradient amplitude constraints
%       are unnecessary.
%

%   ----- Get a list of the amplitude constraints per point -----

if (D==1)           % 1D, so just magnitude constraints.
    if (l1grad==0)
        Ag = [1;-1];
        bg = [gmax;gmax];
    else            % Only need one constraint if type 1 or 2
        Ag=1;        
        bg=gmax;
    end;
elseif (D==2)
    [Ag,bg] = pwlapprox(gmax, npwl,1);
    if (l1grad==1)      % Only keep quadrant 1 constraints
        s=size(Ag);
        Ag = Ag(1:s(1)/4,:);
        bg = bg(1:s(1)/4,:);
    end;
elseif (D==3)
    [Ag,bg] = pwl3d9approx(gmax);
    if(type==0)
        error('Gmax not implemented for type 0 and 3D problems');
    end;
end;
consperpoint = length(bg);


%   ------ Next, some constraints are redundant based on max
%       slew constraints and boundary constraints -------

igstart = 0*bg;     % # points at start where can ignore pth constraint.
igend = 0*bg;       % # points at end...
if (l1grad==1)
    g0test=abs(g0);     % need to put g0 in quadrant 1.
    gftest=abs(gf);
else
    g0test=g0;  
    gftest=gf;
end;


if (1==1)   % Eliminate Redundant Gradient constraints.
        % These are constraints that could not be violated
        % since the boundary point is too far for the constraint to be
        % reached in the given # of points.

    bestsmax = eta*Vmax/min(Lcoil) + max(Rcoil)/min(Lcoil)*gmax;
        % Use highest possible slew rate for test!

    for p=1:consperpoint;
    if (~isempty(g0))
        g0dist = normdist(Ag(p,:).',bg(p,1),g0test.');
        igstart(p) = floor(abs(g0dist/T/bestsmax)); 
    end;
    if (~isempty(gf))
        gfdist = normdist(Ag(p,:).',bg(p,1),gftest.');
        igend(p) = floor(abs(gfdist/T/bestsmax));   
    end;
    end;
end;
igend = N+1-igend;  % Convert to point number.

    
%   ------ Now expand the constraints to apply to each point
%   
ngcons=0;
Agcons=zeros(N*consperpoint,N*D);
bgcons=zeros(N*consperpoint,1);


for k=1:N
    for p=1:consperpoint
    if ((k >= igstart(p)) & (k <= igend(p)))    % be conservative!
        ngcons=ngcons+1;
        for d=1:D
            Agcons(ngcons,(d-1)*N+k)=Ag(p,d);   % be conservative!
        end;
        bgcons(ngcons,1)=bg(p); 
    end;
    end;
end;
Agcons = Agcons(1:ngcons,:);    % Discard zero rows at end.
bgcons = bgcons(1:ngcons,:);    % Discard zero rows at end.



% =============================================================
%   Slew Constraints
% =============================================================
%
%   Note 4: These are similar to the gradient constraints (Note 3).
%       However, there are only N-1 "points", as the DIFFERENCE
%       between gradient points is what is constrained.
%
%       Similar to gradient constraints, these are |g[k]-g[k-1]|
%       for 1D, but quadratic for 2D or 3D.  Again, the quadratic
%       constraints are approximated by linear.  
%
%       If type=2, the slew constraints can also be limited to
%       those that have positive coefficients on g[k].
%       
%       Slew constraints can't be reduced based on boundary 
%       conditions though.
 

        
%   ------ First get a list of the constraints per point -----

if (D==1)
    if (l1slew==0)
        As = [1;-1];
        bs = [1; 1]*eta*Vmax*T/Lcoil;
    else
        As = 1;
        bs = eta*Vmax*T/Lcoil;
    end;
elseif (D==2)
    [As,bs] = pwlapprox(eta*Vmax*T/Lcoil, npwl,1);
    if (l1slew==1)      % Only keep quadrant-1 constraints
        s=size(As);
        As = As(1:s(1)/4,:);
        bs = bs(1:s(1)/4,:);
    end;
elseif (D==3)
    [As,bs] = pwl3d9approx(eta*Vmax*T/Lcoil);
    if (type<2)
        error('Smax constraints and type=0,1 not implemented');
    end;
end;
consperpoint = length(bs);


%   ----- Now expand these for N points. -----

Avscons = zeros((N-1)*consperpoint, N*D);
Avgcons = zeros((N-1)*consperpoint, N*D);
Avmcons = zeros((N-1)*consperpoint, (N-1)*D);
bscons = zeros((N-1)*consperpoint, 1);
RoverL = Rcoil/Lcoil;   
for k=1:(N-1)
    for p=1:consperpoint
    for d=1:D
        Avscons((k-1) *consperpoint+p,(d-1)*N+k  )=    -As(p,d);
        Avscons((k-1) *consperpoint+p,(d-1)*N+k+1)=     As(p,d);
        Avgcons((k-1) *consperpoint+p,(d-1)*N+k  )=     As(p,d)*RoverL*T/2;
        Avgcons((k-1) *consperpoint+p,(d-1)*N+k+1)=     As(p,d)*RoverL*T/2;
        Avmcons((k-1)*consperpoint+p,(d-1)*(N-1)+k )=   As(p,d);
    end;
        bvscons((k-1)*consperpoint+p,1)=bs(p);
    end;
end;



% ============================================================
%       Combine All Constraints
% ============================================================
%
%   Note 5: Here it is worth looking back at Notes 1 and 2.
%   The constraints are all put together to constrain the
%   X variable.



%   ---------- Type 0 ----------

if ((l1grad==0) & (l1slew==0))
    A = [Ag0; Agf; Am; Agcons; Avscons+Avgcons];
    b = [bg0; bgf; bm; bgcons; bvscons];
    s = size(A); 
    f = zeros(s(2),1);
    f(1:N)=1;
end;


%   ---------- Type 1 ----------

if ((l1grad==1) & (l1slew==0))
    
    % --- Form constraints for "h" variables, h>|g| (Note 2) ---
    disp('Type 1 Formulation');

    I = eye(N*D);
    Agm = [I -I ; -I -I ];
    Agm = [eye(N*D) -eye(N*D);-eye(N*D) -eye(N*D)];
    bgm = zeros(N*D*2,1);


    % --- Combine all constraints ---

    A = [Ag0 0*Ag0; Agf 0*Agf; Am 0*Am;];
    A = [A; 0*Agcons Agcons; Avscons+Avgcons 0*Avscons; Agm];
    b = [bg0; bgf; bm; bgcons; bvscons; bgm];
    f = [zeros(N*D,1); ones(N*D,1)];
end;

if ((l1grad==1) & (l1slew==1))

    % --- Form constraints for "h" variables, h>|g| (Note 2) ---
    disp('Type 2 Formulation');

    I = eye(N*D);
    Z = zeros(N*D,(N-1)*D);
    Agm = [I -I Z; -I -I Z];
    bgm = zeros(N*D*2,1);


    % --- Form constraints for "s" variables, 
    %           
    %   s> | dg/dt + R*g/L |    (modified for voltage limit)
    %

    %   (N-1 of these per dimension)

    I1 = (RoverL*T/2-1)*eye(N-1,N) + (RoverL*T/2+1)*[zeros(N-1,1) eye(N-1)]; 
    A2 = eye((N-1)*D);
    if (D==1)
        A1=[I1];
    elseif (D==2)
        A1=[I1 0*I1; 0*I1 I1];
    elseif (D==3)
        A1=[I1 0*I1 0*I1; 0*I1 I1 0*I1; 0*I1 0*I1 I1];
    end;
    Asm = [ A1 0*A1 -A2; -A1 0*A1 -A2];
    bsm = zeros(2*(N-1)*D,1);


    % --- First assemble A as if "s" variables aren't there. ---
    A = [Ag0 0*Ag0; Agf 0*Agf; Am 0*Am;];

    % --- Zero pad on right side, as the above don't apply to "s" vars. ---
    s = size(A);
    A = [A zeros(s(1),(N-1)*D)];

    % --- Add in 1st-quadrant gradient and slew magnitude constraints ----
    s1 = size(Avmcons);
    s2 = size(Agcons);
    A = [A; 0*Agcons Agcons zeros(s2(1),(N-1)*D)];
    A = [A; zeros(s1(1),2*N*D) Avmcons];

    % --- Add in constraints to force |g| and |s| to converge. ---
    A = [A; Agm; Asm];

    % --- b is just the stack of all these... ----
    b = [bg0; bgf; bm; bgcons; bvscons; bgm; bsm];

    % --- f places cost on h and s variables.
    f = [zeros(N*D,1); ones(N*D,1); ones((N-1)*D,1)];

    % ---- Try to get it to end earlier!! ----

        if (1==0)   % polynomial-increasing cost
        fm = ([1:N]/N).^2;
        fh = ones(D,1)*fm;
        fs = ones(D,1)*fm(1:N-1);
        f = [zeros(N*D,1); fh(:); fs(:)];
        end;


        if (1==0)   % step increase in cost.
            te=10;
            tm=1000;
            for dd=1:D
            tf=N*D+N*dd;
            f(tf-te:tf)=tm*f(tf-te:tf);
            tf=2*N*D+(N-1)*dd;
            f(tf-te:tf)=tm*f(tf-te:tf);
            end;
        end;
end;
    

s=size(A);
tt = sprintf('Calling linprog() with %d variables and %d constraints',s(2),s(1));
disp(tt);


if (1==0)
    f
    A
    b
end;


if (exist('linprog')~=0)
  [g,fval,v] = linprog(f,A,b);
  if (v==0)
    lpopt = optimset('MaxIter',500);
    [g,fval,v] = linprog(f,A,b,[],[],[],[],[],lpopt);
  end;
else
  sz = size(A);
  g = zeros(sz(2),1);
  fval =0;
  v = 1;
  disp('*** linprog does not exist.  Returning zero solution. ***');
end;

g = reshape(g(1:N*D),N,D);


% --------- Record End time, and display times ---------

n2 = clock; 
if (0==1)
    tt=sprintf('Start Time  %02d:%02d:%02d ',round(n1(4:6)));
    disp(tt);
    tt=sprintf('End Time    %02d:%02d:%02d ',round(n2(4:6)));
    disp(tt);
end;
