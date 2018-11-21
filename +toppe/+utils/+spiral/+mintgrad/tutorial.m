%
%	Full tutorial for minimum-time gradient waveform
%	design website.
%


clear;
close all;
ex1=1;
ex2=1;
ex3=1;

disp('-----------------------------------------------------------------------');
disp('Minimum-Time Multi-Dimensional Waveform Desing with Convex Optimization');
disp('-----------------------------------------------------------------------');
disp(' ');
disp(' ');
if (ex1)	% Flag to switch off.
disp('EXAMPLE 1. ');
disp('---------- ');
disp('We begin by looking at the design of a simple gradient lobe, with');
disp('constant limits on gradient amplitude and slew rate.  ');
disp('It is relatively obvious that the solution to this problem is ');
disp('a trapezoidal waveform, but this should help us understand what ');
disp('the function is doing. ');
disp(' ');
disp('Let us assume that we simply want a 1D gradient that has area ');
disp('of 8 cm^(-1) in k-space (corresponding to .6mm resolution),  ');
disp('which starts and ends at zero-amplitude.  ');
disp('socpgrad.m will tell us if there is a solution of length N:' );
disp(' ');
disp('T=0.00002;	% Gradient Sampling period of 20us. ');
T=0.00002;	
disp('N=20;');
N=20;
disp('[g,v] = socpgrad(N,0,0,8,T,4,15000);');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');
[g,v] = socpgrad(N,0,0,8,T,4,15000);
disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ');
disp('Okay, there was no solution for N=20 (ie 400us duration).');
disp(' ');
disp('Now we will try with N=50: ');
disp(' ');
disp('N=50;');
N=50;
disp('[g,v] = socpgrad(N,0,0,8,T,4,15000);');
disp('<<Press Any Key>>'); pause;
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');
[g,v] = socpgrad(N,0,0,8,T,4,15000);
disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ');
disp('This time there is a solution.  We will plot it:');
disp('plotgradinfo(g,T);');
disp('<<Press Any Key>>'); pause; disp(' ');
plotgradinfo(g,T);

disp(' ');
disp('plotgradinfo() plots some information about the gradient waveforms. ');
disp('It is important to pass the sampling rate to this function! ');
disp('Notice (top center and top right) that all the desired constraints ');
disp('are met:  the gradient starts and ends at 0, and has area 5cm^(-1). ');
disp(' ');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('However, the gradient does not look trapezoidal.  We still need ');
disp('to check if a solution exists for lower N.  We can do this until ');
disp('N is too small to return a solution: ');
disp(' ');
disp('v=0;');
disp('N=50;');
disp('while (v>0)');
disp('[g,v] = socpgrad(N,0,0,8,T,4,15000);');
disp('if v>0 grad=g; N=N-1; end;');
disp('end;');
disp('plotgradinfo(g,T);');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');

N=50;
v=1;
while (v>0)
  [g,v] = socpgrad(N,0,0,8,T,4,15000);
  if v>0 grad=g; N=N-1; end;
end;
plotgradinfo(grad,T);
disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');

disp(' ');
disp('So a solution exists with N=38, but not with N=37.  ');
disp('Since the convex optimization functions are guaranteed to converge,');
disp('we can assume that N=38 is the minimum-time solution given the');
disp('constraints.  Notice that the waveform (top right) does look');
disp('trapezoidal now. ');

disp(' ');
disp('==================================================================== ');
disp('<<Press Any Key>>'); pause; disp(' ');
disp(' ');
end;
if (ex2)	% Flag to switch off.
disp('EXAMPLE 2. ');
disp('----------');
disp(' ');
disp('Now we begin to look at 2D gradient design.  Here we will ');
disp('assume that we want the gradients to be freely rotatable  ');
disp('for oblique scanning.  Thus the amplitude and slew rate ');
disp('limits are quadratic. ');
disp(' ');
disp('First, we will now use mintimegrad(), instead of just ')
disp('searching solutions until N is small enough.  mintimegrad() finds');
disp('the optimal solution length much more quickly, by bisecting intervals');
disp('that must contain the optimal length.');
disp(' ');
disp('We will repeat the design of Example 1, but with  ');
disp('a desired k-space change of 7cm^(-1) in x and 6cm^(-1) ');
disp('in y: ');
disp('[g] = mintimegrad(50,[0 0],[0 0],[7 6],T,4,15000,0,3);');
disp('plotgradinfo(g,T);');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');

T=0.00002;	
[g] = mintimegrad(50,[0 0],[0 0],[7 6],T,4,15000,0,3);
plotgradinfo(g,T);

disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ');
disp('So now the solution is of length 42, or 840us.')
disp('Not surprisingly, the solution is still a trapezoid,');
disp('and the gradient/slew rate is simply allocated between');
disp('x and y to satisfy the maximum constraint.');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('However, if we stipulate that the x-gradient must ');
disp('start with a value of 4 G/cm, the solution is no longer');
disp('obvious.  This is where the power of the convex optimization');
disp('method starts to show:');
disp(' ');
disp('[g] = mintimegrad(50,[4 0],[0 0],[7 6],T,4,15000,0,3);');
disp('plotgradinfo(g,T);');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');

  [g] = mintimegrad(50,[4 0],[0 0],[7 6],T,4,15000,0,3);
  figure(2);
  plotgradinfo(g,T);

disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ');
disp('Compare the solutions in Figure 1 and Figure 2.   ');
disp(' ');
disp('Notice that the waveform in Figure 2 is considerably ');
disp('more complex.  The waveform generally reaches either ');
disp('the amplitude or slew-rate limit over its entire ');
disp('duration, as you would hope. ');
disp(' ');
disp('<<Press Any Key>>'); pause; disp(' ');
disp(' ');
disp('As a final example here, what if we ALSO want the  ');
disp('gradient to be flow-compensated?  That is, we want ');
disp('the first moments to be zero at the end of the  ');
disp('waveform?  This can easily be added into the ');
disp('constraints of the waveform design: ');
disp(' ');
disp('[g] = mintimegrad(60,[4 0],[0 0],[7 6;0 0],T,4,15000,0,3);');
disp('plotgradinfo(g,T);');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');
  [g] = mintimegrad(60,[4 0],[0 0],[7 6;0 0],T,4,15000,0,3);
  plotgradinfo(g,T);
disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ');
disp('Notice (Figure 1) that the first moment (bottom center) now ');
disp('ends at zero, and that the solution length is considerably ');
disp('longer - 1.5ms instead of 720us.');
disp(' ');
disp(' ');
disp(' ');
disp('==================================================================== ');
disp('<<Press Any Key>>'); pause; disp(' ');
disp(' ');
end;
if (ex3)	% Flag to switch off.
disp('EXAMPLE 3. ');
disp('----------');
disp(' ');
disp('This example shows an important application of optimal-time ');
disp('gradient design.  Here we calculate the rewinder gradients ');
disp('for spiral imaging trajectories. ');
disp(' ');
disp('First, we design the spiral trajectory ');
disp('(with FOV=16cm, 1mm resolution, 60interleaves): ');
disp(' ');
disp('[k,g] = vds(15000,4,T,60,16,0,0,5) ');
disp('<<Press Any Key>>'); pause; disp(' ');
disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');
close all;
T=.00002;
[k,g] = vds(15000,4,T/5,60,16,0,0,5);
g = g(1:5:length(g));
disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');
disp(' ')
disp('Now we calculate the area and moments to rewind, and call ');
disp('mintimegrad() again: ')
disp('gsp = [real(g(:)) imag(g(:))]; ');
disp('[g,k,s,m1] = calcgradinfo(gsp,T); ');
disp('sz=size(k); ');
disp(' ')
disp('g=mintimegrad(60,gsp(sz(1),:),[0 0],[-k(sz(1),:);-m1(sz(1),:)],T,4,15000,0,3);');
disp(' ')
disp('<<Press Any Key>>'); pause; disp(' ');
  disp('- - - - - - - - - (Start of Function Output)  - - - - - - - - - - - ');
  gsp = [real(g(:)) imag(g(:))];
  [k,g,s,m1] = calcgradinfo(gsp,T);
  sz=size(k);
  t0=T*sz(1);
g=mintimegrad(60,gsp(sz(1),:),[0 0],[-k(sz(1),:);-m1(sz(1),:)],T,4,15000,t0,3);
  plotgradinfo([gsp;g],T);
  disp('- - - - - - - - - (End of Function Output)  - - - - - - - - - - - ');

disp(' ');
disp('Here we plotted the spiral gradient as well as the rewinder.');
disp('We could have done the rewinder without rewinding the first-moment');
disp('by leaving out the -m1(sz(1),:) term.');
end;



