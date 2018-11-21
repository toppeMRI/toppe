function [x,info,z,w,hist,time]= ...
	socp(f,A,b,C,d,N,x,z,w,abs_tol,rel_tol,target,max_iter,Nu,out_mode)

% Solve Second-order Cone Problem by primal-dual interior point method.
%
%   [x,info,z,w,hist,time] = ...
%       socp(f,A,b,C,d,N,x,z,w,abs_tol,rel_tol,target,max_iter,Nu,out_mode)
%
%   Input and/or output parameters may be omitted, starting from the end.
%   For input parameters, default values are then used (or, in the cases
%   of x and of z,w, a procedure is used to compute them).  This is also
%   done when a parameter is the empty list, [].  The shortest calling
%   sequence is:
%
%   x = socp(f,A,b,C,d,N)
%
% PROBLEM:
%   minimize      f'*x
%   subject to    norm(A_i*x+b_i) <= c_i'*x+d_i ,  i=1,...,L
%
%       The dual problem is:
%           maximize      -(b'*z+d*w)
%           subject to    A'*z+C'*w == f
%                         norm(z_i) <= w_i ,  i=1,...,L
%       (see description of input arguments for definition of A,b,C,d,z,w)
%
% INPUT ARGUMENTS:
%
%  PROBLEM DESCRIPTION:
%   f -- vector defining linear objective, length(f) == length(x)
%   A, b, C, d -- constraints:
%       A is a matrix with the A(i) vertically stacked:
%          A = [ A(1); A(2); ... A(L)]
%       b is a vector with the b(i) vertically stacked:
%          b = [ b(1); b(2); ... b(L)]
%       C is a matrix with the c(i)' vertically stacked:
%          A = [ c(1)'; c(2)'; ... c(L)']
%       d is a vector with the d(i) vertically stacked:
%          d = [ d(1); d(2); ... d(L)]
%       The matrix [A;C] should be full rank; any problem can be
%       reformulated to satisfy this condition.
%   N -- vector of size L, defining the size of each constraint,
%       i.e.  size(A(i)) == [N(i),n]
%             size(b(i)) == [N(i),1]
%       For linear constraints define N(i)=0 (i.e. A(i)=[], b(i)=[])
%           and the constraint will be:  0 <= c(i)'*x+d(i)
%
%  INITIAL POINT:
%   x -- strictly primal feasible initial point. Must satisfy:
%            norm(A(i)*x+b(i)) < c(i)'*x+d(i) ,  i=1,...,L
%   z,w -- strictly dual feasible initial point. We define:
%            z = [z(1); ... ; z(L)]
%            w = [w(1); ... ; w(L)]
%       where  size(z(i))==[N(i),1] and w(i) are scalars.
%       The dual feasibility conditions are:
%            A'*z + C'*w == f
%            norm(z(i)) < w(i) ,  i=1,...,L
%       If either x or z,w are omitted or empty (i.e. x=[], or z=[] and w=[]),
%       an appropriate method is used to find a strictly feasible initial
%       point (phase 1 / big-M)
%
%  STOPPING CRITERIA: (stop when any of the following is met)
%   abs_tol -- maximum absolute error in objective function;
%       guarantees that for any x:  abs(f'*x - f'*x_opt) <= abs_tol
%   rel_tol -- maximum relative error in objective function;
%       guarantees that for any x:
%       abs(f'*x - f'*x_opt)/(f'*x_opt) <= rel_tol  (if f'*x_opt > 0)
%       Negative value has special meaning, see target
%   target -- if rel_tol<0, stops when f'*x<target or -b'*z>=target
%   max_iter -- limit on number of algorithm outer iterations.
%       Most problems can be solved in less than 50 iterations.
%       Called with max_iter=0 only checks feasibility of x and z,
%	(and returns gap and deviation from centrality).
%
%  OTHER PARAMETERS:
%   Nu -- duality gap vs. deviation from centrality reduction weight.
%       As a general rule, larger values of Nu yield faster convergence,
%       up to a point where the deviation from centrality becomes too
%       large and the convergence becomes very slow.
%       Required Nu>0, recommended Nu>1, typical range 2<=Nu<=50. Try Nu=10.
%   out_mode -- specifies what will be output in hist:
% 	0: hist is returned empty (default value)
%	1: vector with duality gap (an upper bound on absolute error),
%	   for the initial point and after each iteration
%	2: matrix with duality gap and deviation from centrality
%
% OUTPUT ARGUMENTS:
%   x -- solution
%   info -- string:
%           'absolute accuracy reached', 'relative accuracy reached',
%           'target value reached', 'target value is unachievable',
%           'maximum iterations exceeded', 'error'
%   z, w -- solution to the dual problem
%   hist -- see out_mode
%   time -- statistics, vector with 3 numbers:
%       time(1) = user time (cpu time used)
%       time(2) = system time (time spent in system calls)
%       time(3) = total number of iterations performed
%

% 1997, mlobo@isl.stanford.edu


%%% add equality constraints?...


%%% CONSTANTS

check_rank=0;	% 1 to do rank reduction procedure, 0 to skip it
check_feas=0;	% 1 to do check of initial point feasibility, 0 to skip it

eq_eps_in=1e-10;	% to decide when initial point is not in hyperplane
eq_eps_out=1e-4;	% same for final point
			%  (these should probably be made to depend on m)

BigM_K=2;	% scaling factor for expanding big-M bound
BigM_iter=2;	% number of iterations before deciding whether to expand


%%% DEFAULT PARAMETERS

Nin=6;	% min. num. of parameters

if nargin<Nin+9,		out_mode=[];
 if nargin<Nin+8,		Nu=[];
  if nargin<Nin+7,		max_iter=[];
   if nargin<Nin+6,		target=[];
    if nargin<Nin+5,		rel_tol=[];
     if nargin<Nin+4,		abs_tol=[];
      if nargin<Nin+3,		w=[];
       if nargin<Nin+2,		z=[];
        if nargin<Nin+1,	x=[];
         if nargin<Nin,	error('insuficient number of parameters');
end; end; end; end; end; end; end; end; end; end;

if isempty(abs_tol),	abs_tol=1e-6;	end;
if isempty(rel_tol),	rel_tol=1e-4;	end;
if isempty(target),	target=0;	end;
if isempty(max_iter),	max_iter=100;	end;
if isempty(Nu),		Nu=10;		end;
if isempty(out_mode),	out_mode=0;	end;


		% socp makes a recursive call for phase 1
		% with the data structure already changed;
		% this is signaled by an empty d;
if isempty(d),

		% the following flag tells socp not to change back
		% to the entry structure before exiting
	abcd=0;

	%%% SIZES

	L=length(N);
	[m,n]=size(A);

else

	%%% SIZES

	% make all vectors column vectors
	f=f(:); b=b(:); d=d(:); N=N(:); x=x(:); z=z(:); w=w(:);

	n=length(f);
	if ~(sum(N)==0) & ~(size(A,2)==n),
		error('Number of columns of A must equal length of f');  end;
	if ~(size(C,2)==n),
		error('Number of columns of C must equal length of f');  end;
	if ~isempty(x) & ~(length(x)==n),
		error('f and x must have the same length');  end;

	L=size(C,1);
	if ~(length(N)==L),
		error('length of N must equal number of rows in C');  end;
	if ~(length(d)==L),
		error('length of d must equal number of rows in C');  end;
	if ~isempty(w) & ~(length(w)==L),
		error('length of w must equal number of rows in C');  end;

	m=size(A,1);
	if ~(sum(N)==m),
		error('sum(N) must equal number of rows in A');  end;
	if ~(length(b)==m),
		error('sum(N) must equal length of b');  end;
	if ~isempty(w) & ~(length(z)==m),
		error('sum(N) must equal length of z');  end;
	m=m+L;


        %%% CHANGE TO THE DATA STRUCTURE USED IN THE C-FCT

	abcd=1;
	N=N+1;
		% build permutation matrix
	i=(1:m)';
	j=[];
	for k=1:L,
		j=[j; (sum(N(1:k-1))+1:sum(N(1:k))-1)'];
	end;
	j=[j; cumsum(N)];
	P=sparse(i,j,ones(m,1),m,m,m);
		% permute
	A=P'*[A;C];
	b=P'*[b;d];
	if isempty(w),
		z=[];
	else
		z=P'*[z;w];
	end
end;


%%% BUILD CONE VARIABLE SELECTORS

SU=sparse([],[],[],m,L,m-L);
ST=sparse([],[],[],m,L,L);
T=zeros(m,1);
k=1;
for i=1:L,
 if N(i)>1, SU(k:k+N(i)-2,i)=ones(N(i)-1,1); end;
 ST(k+N(i)-1,i)=1;
 T(k+N(i)-1)=1;
 k=k+N(i);
end;


%%% RANK REDUCTION

rank_red=0;
if check_rank,
    [U,S,V]=svd(A,0);
    r=rank(S);
    if r<n,
	t=norm(V(:,r+1:n)'*f);
	if t>10*eps,
	    disp(['Warning: if feasible, original problem is unbounded; ',...
		'f is not orthogonal to null-space of A; ',...
		'norm(V2''*f)=', num2str(t)]);
	end;
	V1=V(:,1:r);
	frr=f;
	f=V1'*f;
	A=U(:,1:r)*S(1:r,1:r);
	x=V1'*x;
	disp(['Columns of A are not independent, ',...
		int2str(n-r),' variable(s) removed.']);
	rank_red=1;
    end;
end;


%%% PHASE 1

if isempty(x),

     % To define minimum initial slack
     % we take a guess at the order of mag. of u
     % (a bad guess makes for a very off-center initial point
     % so this may need to be edited for specific problems)
    slacka=0.1*norm([A b], inf);
     % PICK ANY PRIMAL POINT
    x=zeros(n,1);
    u0=A*x+b;				% point in feasible hyperplane
    k=-ST'*u0+sqrt(SU'*(u0.^2));	% cone "infeasibility" vector
    alpha=slacka+max(k);		% to make all strictly feas.
    x0=A\T;
    if norm(T-A*x0)<10*eps*sqrt(m),	% XXX factor of 10 is arbitrary
         % A FEASIBLE POINT CAN BE TRIVIALLY FOUND
        x=alpha*x0;
    else
	if alpha>=0,			% if alpha<0, skip phase 1
	     % EXTENDED PROBLEM
	    x=[x; alpha];		% extend primal variable
	    f1=[zeros(n,1); 1];		% phase 1 cost
	    A1=[A T];			% extend u: t(i)+=alpha
	     % SOLVE FEASIBILITY PROBLEM
	     % i.e. stop when alpha<0 (or dual obj.>=0)
	    [x,info1,z1,dummy,hist1,time1]= ...
		    socp(f1,A1,b,[],[],N,x,[],[],0,-1,0,max_iter,Nu);
	    if x(n+1)>=0,
		    error(['Phase 1 failed, alpha>=0;  info=',info1]);
	    end;
	    x=x(1:n);		% remove extra variable from solution
	end;
    end;

elseif check_feas,	% CHECK FEASIBILITY OF PRIMAL INITIAL POINT

    u=A*x+b;
    k=0;
    for i=1:length(N),
	if N(i)>1,
	    if norm(u(1+k:N(i)-1+k))>=u(N(i)+k),
		    error(['Primal point is not strictly feasible: ',...
		     'constraint ',int2str(i),' is not satisfied.']);
	    end;
	elseif N(i)==1,
	    if u(1+k)<=0,
		    error(['Primal point is not strictly feasible: ',...
		     'constraint ',int2str(i),' is not satisfied.']);
	    end;
	else
	    error(['Invalid dimension: N(',int2str(i),')<1']);
	end;
	k=k+N(i);
    end;

end;


%%% PHASE 2  (also phase 1 in recursive call)

if isempty(z),

    %%% BIG-M (with expanding bound)

    BigM_repeat=max_iter/BigM_iter;
     % To define minimum initial slack
     % we take a guess at the order of mag. of z
     % (a bad guess makes for a very off-center initial point
     % so this may need to be edited for specific problems)
    slackb=0.1*norm(f,inf)/norm(A, 1);
     % PRIMAL BOUND
    u=A*x+b;			% point in feasible hyperplane
    B0=BigM_K*(T'*u);		% bound on sum(t)
     % DUAL FEASIBLE POINT
    z0=A'\f;			% pick any point satisfying hyperplane constr.
    k=-ST'*z0+sqrt(SU'*(z0.^2));	% cone "infeasibility" vector
    beta=slackb+max(0,max(k));	% to make all strictly feas. (incl. beta>0)
    z=z0+T*beta;			% add beta to all w(i)
    z1=[z; beta];			% extend dual variable
     % EXTENDED PROBLEM
    A1=[A; -T'*A];			% extend u:  tB=B0-sum(t)
    b1=[b; (B0-T'*b)];
    N1=[N(:); 1];
     % SOLVE PROBLEM
    if exist('time1'),
    	hist=hist1;
        time=time1;
    else
        hist=[];
        time=[0,0,0];
    end;
    i=0;
    info=0;
    while info==0,
	i=i+1;
	[x,info,z1,hist1,time1]= socp_mex(f,A1,b1,N1,x,z1, ...
		abs_tol,rel_tol,target,BigM_iter,Nu,out_mode);
	time=time+time1;
	hist=[hist, hist1(:,(i>1)+1:time1(3)+1)];
	if info==4 & i<BigM_repeat,
	    info=0;
	end;
	u=A*x+b;
	t=T'*u;
	if BigM_K*t>B0,
	    b1(m+1)=b1(m+1)+BigM_K*t-B0;
	    B0=BigM_K*t;
	end;
    end;
     % check quality of dual solution
     % (choice of 1e-4 for displaying warnings is arbitrary)
    e1=z1(m+1)/norm(z1(1:m));
    if rel_tol>=0 & e1>max(rel_tol,1e-4),
     disp('Warning: big M bound appears to be active, solution may be wrong.');
     disp([' Relative dual slack from bound: ',num2str(e1)]);
    end;
    z=z1(1:m);	% remove extra variable from solution

else

    %%% SOLVE WITHOUT BIG-M  (initial z given)

    if check_feas,		% check feasibility of dual initial point
	e1=norm(A'*z-f)/norm(f);
	if e1 > eq_eps_in,
	    error(['Dual point is not feasible.  ',...
	     'Rel. error in equality constraint=',...
	     num2str(e1)]);
	end;
	k=0;
	for i=1:length(N),
	    if N(i)>1,
		if norm(z(1+k:N(i)-1+k))>=z(N(i)+k),
			error(['Dual point is not strictly feasible: ',...
			 'constraint ',int2str(i),' is not satisfied.']);
		end;
	    elseif N(i)==1,
		if z(1+k)<=0,
			error(['Dual point is not strictly feasible: ',...
			 'constraint ',int2str(i),' is not satisfied.']);
		end;
	    else
		error(['Invalid dimension: N(',int2str(i),')<1']);
	    end;
	    k=k+N(i);
	end;
    end;

    [x,info,z,hist,time]= socp_mex(f,A,b,N,x,z, ...
		abs_tol,rel_tol,target,max_iter,Nu,out_mode);

end;


%%% UNDO RANK REDUCTION TRANSFORMATION

if rank_red,
    f=frr;
    x=V1*x;
end;


%%% DIAGNOSTICS

if info==1 | info==2,
    e2=norm(A'*z-f)/norm(z);
    if e2 > eq_eps_out,
	disp('Dual is not in the feasible hyperplane, solution may be wrong.');
	disp(['  Relative error in equality contraint=',num2str(e2)]);
    end;
end;


%%% CONVERT INFO TO STRING

if	info==1,	info='absolute accuracy reached';
elseif	info==2,	info='relative accuracy reached';
elseif	info==3,
	if f'*x<=target,
			info='target value reached';
	else		info='target value is unachievable';
	end;
elseif	info==4,	info='maximum iterations exceeded';
else			info='error';
end;


%%% CHANGE z BACK TO z,w

if abcd,
	z=P*z;
	w=z(m-L+1:m);
	z=z(1:m-L);
else
	w=[];
end;
