function [g,k,t,s,npts] = genspiralvd(D, N, nl, gamp, gslew, dsamp) 
% /*
%  *
%  *    Subroutine for generating slew-rate-limited spirals.  
%  *    Returns length of resulting gradients, after adding downramp.
%  *
%  *    Craig Meyer, Copyright Leland Stanford Junior University, 1996.
%  *
%  *    variable density added by D. Noll, 10/25/01
%  */
% 
% #define MAXDECRATIO 32 /* maximum allowed decimation of input ts */
% #define ABS(a) ( (a > 0) ? a : -a ) 
MAXDECRATIO=32;
%     double A, risetime,GAM,OM,S,absk,targetk;
%     short delx, dely;
%     int m, n, npts, dnpts, loop, res;
%     int dentrans, den1, kdenrad;
%     double decratio, om, s;
%     double g0, thetan_1, theta, deltheta, taun_1, taun, tauhat;
%     double absg,gtilde,B,t1,t2,t3,ac,tgx,tgy;
%     double kx, ky, oldkx, oldky;
%     double OMF, omf, denrad, scoffset, denoffset, scthat,fractrans, realn, ksv;


% npts is the number of samples in the spiral itself and not the shootout
% grad.
%%%%%%%%%%%%%% Initialize scanner variables  %%%%%%%%%%%%%%%%%%%%%
GRESMAX= 21000;
if ~exist('nl','var')
    nl=1   % factor of undersampling in the outer regions of kspace
end
if ~exist('gamp','var')
      gamp=2.2; %3.50; % 2.2 for both 1.5 T and 3 T data
end
if ~exist('gslew','var')
      gslew=180 % 200 % 180 for 3T data and 120 (150) for 1.5 T data
end
nramp=0;
% nramp=100;
MAX_PG_WAMP=32766;

gts = 4e-06;

Tmax = GRESMAX*gts;

dts = gts;
opfov = D;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GAM = 4257.0;
A = MAX_PG_WAMP;
risetime = gamp/gslew*10000;
OM = 2.0*pi/nl*D/(1/(GAM*gamp*dts));
S = (dts/1e-6)*A/risetime;
absk = 0.;;
targetk = N/2;
OMF = OM*nl;
dentrans = dsamp/2;
S0=gslew*100;


ac = A;

loop = 1;
decratio = 1;

while (loop)
    loop = 0;
    om = OM/decratio;
    s = S/decratio;
    omf = OMF/decratio;
    den1 = 0;
    g0 = 0;
    Gx(1) = g0;
    Gy(1) = 0;
    absg = abs(g0);      % absg = hypot(g0,0);
    oldkx = 0;
    oldky = 0;
    kx = Gx(1);
    ky = Gy(1);
    thetan_1 = 0;
    taun = 0;
    n = 1;
    while (absk < targetk)  %keep generating waveform until reach desired k-space distance
        taun_1 = taun;
        taun = sqrt(kx.^2+ky.^2)/A;  %taun = hypot(kx,ky)/A;
        tauhat = taun;
        realn = (n-1)/decratio; %realn = n/decratio
        if (realn >  dsamp)  %undersampled region
            if (den1 == 0)  
                ksv = taun;
                den1 = 1;
            end
            if (realn >  (dsamp+dentrans)) 
                scoffset = scthat;
                denoffset = taun_1;
                scthat = scoffset + om*(tauhat - denoffset);
                fractrans = 1;
            else 
                scoffset = scthat;
                denoffset = taun_1;
                fractrans = (realn - dsamp)/( dentrans); 
                fractrans = 1 - ( (fractrans-1)*(fractrans-1));
                scthat = scoffset + (omf + (om-omf)*fractrans)*(tauhat - denoffset);
            end
        else   %fully sampled region
            scthat = omf*tauhat;
            fractrans = 0;
        end
        
        theta = atan2(scthat,1.0)+scthat;
        if (absg < ac)
            deltheta = theta-thetan_1;
            B = 1.0/(1.0+tan(deltheta)*tan(deltheta));
            gtilde = absg;
            t1 = s*s;
            t2 = gtilde*gtilde*(1-B);
            if (t2 > t1)
                decratio = decratio * 2.0;
                if (decratio > MAXDECRATIO)
                    error('Genspiral failed');
                    return;
                end
                loop = 1;
                break;
            end
            t3 = sqrt(t1-t2);
            absg = sqrt(B)*gtilde+t3;
            if (absg > ac)
                absg = ac;
            end
        end
        
        tgx = absg*cos(theta);
        tgy = absg*sin(theta);
        kx = kx + tgx;
        ky = ky + tgy;
        thetan_1=theta;
        if ~mod(n,round(decratio))  %IRINT(decratio)
            m = n/round(decratio);
            absk = nl*om*taun/(2*pi);
            if (absk > targetk)
                break;
            end
            Gx(m+1) = round((kx-oldkx)/decratio); %&0xfffe;
            Gy(m+1) = round((ky-oldky)/decratio); %&0xfffe;
            oldkx = kx;
            oldky = ky;
        end
        n=n+1;
    end
end

gres1 = m;
npts = m+1;

g = (Gx + i.*Gy)./(MAX_PG_WAMP/gamp);   %slew rate vector
s = diff(g)./(gts*100);  % grad vector
Kx = cumsum([0 Gx])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
Ky = cumsum([0 Gy])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
k = Kx + i.*Ky;  %kspace vector
t = [0:(length(g)-1)]*gts;
matrix = max(abs(k))*2;
maxg = max(abs(g));
maxs = max(abs(s));

kxvec = real(k);
kyvec = imag(k);
sx = real(s);
sy = imag(s);


% return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                       %
%               E N D  O F  F I L E!!!!!!!              %
%                                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

res = ceil(2*nl*om*taun/(2*pi));
kdenrad = ceil(2*nl*om*ksv/(2*pi));
sprintf('resolution = %d, high dens region = %d\n', res,kdenrad);

%/* add downramp */
t1 = atan2(Gy(npts-1), Gx(npts-1))+pi;
delx = round(S*cos(t1));%rint(S*cos(t1));
dely = round(S*sin(t1));%rint(S*sin(t1));
if (delx < 0)
    delx = (delx+1); % & 0xfffe; 
else
    delx = delx; % & 0xfffe;
end
if (dely < 0)
    dely = (dely+1);% & 0xfffe; 
else
    dely = dely; % & 0xfffe;
end
m = abs(Gx(npts-1)/delx);  %ramp length needed for x
n = abs(Gy(npts-1)/dely);  %ramp length needed for y
if(m>n)
    dnpts = npts+n;
else
    dnpts = npts+m;
end
tgx = Gx(npts-1);
tgy = Gy(npts-1);
%add the ramp points to the end of the gradients
for n=npts:dnpts
    tgx = tgx + delx;
    tgy = tgy + dely;
    Gx(n) = tgx;
    Gy(n) = tgy;
end
Gx(n) = 0;
Gy(n) = 0;

%Check if last point is zero.  If not make it zero or ramp down to zero of
%not already very close.  (hmmm, seems redundant... you guys)
% while ((Gx(n) ~= 0)|(Gy(n) ~= 0))
%     
%     if (abs(Gx(n)) < abs(delx))
%         Gx(n+1) = 0;
%     else
%         Gx(n+1) = Gx(n)+delx;
%     end
%     if (abs(Gy(n)) < abs(dely))
%         Gy(n+1) = 0;
%     else
%         Gy(n+1) = Gy(n)+dely;
%     end
%     n=n+1;
% end
gres = n+1;




totx=0;
toty=0;
for j=1:gres-1
    totx=totx+Gx(j)/(MAX_PG_WAMP/gamp);
    toty=toty+Gy(j)/(MAX_PG_WAMP/gamp);
end


areax = sum(Gx(1:gres-1));
areay = sum(Gy(1:gres-1));

%***************NEW REPHASER *********************/
nramp = ceil(risetime/4);
area = sqrt(areax*areax + areay*areay);
maxtri = A*nramp;

if (area > maxtri) %{
  ntop = ceil((area - maxtri)/gamp);
  bx = (-areax/(nramp + ntop));
  by = (-areay/(nramp + ntop));
  %for (j=0; j<nramp; j++)  {
  for j=0:(nramp-1)
    c = j/nramp;
    Gx(n) = bx*c;
    Gy(n) = by*c;
    n = n+1;
  end %}
  %for (j=0; j<nramp; j++)  {
  for j=0:(ntop-1)
    Gx(n) = bx;
    Gy(n) = by;
    n = n+1;
  end %}
  %for (j=0; j<nramp; j++)  {
  for j=0:(nramp-1)
    c = 1 - j/nramp;
    Gx(n) = bx*c;
    Gy(n) = by*c;
    n = n+1;
  end %}
else %} else {
  ntop = 0;
  bx = -areax/(nramp + ntop);
  by = -areay/(nramp + ntop);
  % for (j=0; j<nramp; j++)  {
  for j=0:(nramp-1)

    c = j/nramp;
    Gx(n) = bx*c;
    Gy(n) = by*c;
    n = n+1;
  end %}
  %for (j=0; j<nramp; j++)  {
  for j=0:(nramp-1)

    c = 1 - j/nramp;
    Gx(n) = bx*c;
    Gy(n) = by*c;
    n = n+1;
  end %}

end %}

Gx(n) = 0;
Gy(n) = 0;
gres = n+1;



    Gxr=Gx(end:-1:1);  
    Gyr=Gy(end:-1:1);  
    
    
g = (Gx + i.*Gy)./(MAX_PG_WAMP/gamp);   %slew rate vector
s = diff(g)./(gts*100);  % grad vector
Kx = cumsum([0 Gx])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
Ky = cumsum([0 Gy])*gts*opfov*GAM./(MAX_PG_WAMP/gamp);
k = Kx + i.*Ky;  %kspace vector
t = [0:(length(g)-1)]*gts;
matrix = max(abs(k))*2;
maxg = max(abs(g));
maxs = max(abs(s));

kxvec = real(k);
kyvec = imag(k);
sx = real(s);
sy = imag(s);
    
%[Gx,Gy,Gxr,Gyr,kx,ky]=genspiralvd_rev(22,64,2,2.2,180,300); figure, subplot(121),plot(Gx/MAX_PG_WAMP*gamp),subplot(122),plot(diff(Gx/MAX_PG_WAMP*gamp/10*1000*25));